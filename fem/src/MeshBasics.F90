!*****************************************************************************/
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

MODULE MeshBasics

    USE ElementDescription
    USE BandwidthOptimize
    USE Interpolation
    USE ParallelUtils
    USE Lists
    USe ListMatrix
    USE MeshAllocations
    USE MeshIO
    USE ElementUtils, ONLY : mGetBoundaryIndexesFromParent, &
        NormalDirection, CreateMatrix, TangentDirections, &
        FreeMatrix, Find_Face, Find_Edge
    USE MortarUtils, ONLY : MarkHaloNodes, GeneratePeriodicProjectors
    USE GeometryFitting, ONLY : CylinderFit, SphereFit, TorusFit
    IMPLICIT NONE

CONTAINS






!------------------------------------------------------------------------------
!> This subroutine is used to fill Def_Dofs array of the solver structure.
!> Note that this subroutine makes no attempt to figure out the index of
!> the body, so all bodies are assigned with the same element definition.
!> A similar array of reduced dimension is also filled so as to figure out
!> the maximal-complexity definition over all solvers which use the same
!> global mesh.
!------------------------------------------------------------------------------
   SUBROUTINE GetDefs(ElementDef, Solver_Def_Dofs, Def_Dofs, &
       Def_Dofs_Update, DG)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*), INTENT(IN) :: ElementDef     !< an element definition string
     INTEGER, INTENT(OUT) :: Solver_Def_Dofs(:,:,:) !< Def_Dofs of the solver structure
     INTEGER, INTENT(INOUT) :: Def_Dofs(:,:)        !< holds the maximal-complexity definition on global mesh
     LOGICAL, INTENT(IN) :: Def_Dofs_Update         !< is .TRUE. when the definition refers to the global mesh
     LOGICAL, INTENT(IN) :: DG
!------------------------------------------------------------------------------
     INTEGER, POINTER :: ind(:)
     INTEGER, TARGET :: Family(10)
     INTEGER :: i,j,l,n

     Family = [1,2,3,4,5,6,7,8,9,10]

     ! The default assumption is that the given element definition is applied 
     ! to all basic element families (note that the element sets 9 and 10 are
     ! not included since the explicit choice of the target family is 
     ! a part of the element definition string when the target index is
     ! deduced to be 9 or 10).
     !
     ind => Family(1:8)
     !
     ! If the element family is specified, change the target family 
     !
     IF (SEQL(ElementDef, 'point') )     ind => Family(1:1)
     IF (SEQL(ElementDef, 'line') )      ind => Family(2:2)
     IF (SEQL(ElementDef, 'tri') )       ind => Family(3:3)
     IF (SEQL(ElementDef, 'quad') )      ind => Family(4:4)
     IF (SEQL(ElementDef, 'tetra') )     ind => Family(5:5)
     IF (SEQL(ElementDef, 'pyramid') )   ind => Family(6:6)
     IF (SEQL(ElementDef, 'prism') )     ind => Family(7:7)
     IF (SEQL(ElementDef, 'brick') )     ind => Family(8:8)
     IF (SEQL(ElementDef, 'tri_face') )  ind => Family(9:9)
     IF (SEQL(ElementDef, 'quad_face') ) ind => Family(10:10)

     n = INDEX(ElementDef,'-')
     IF (n<=0) n=LEN_TRIM(ElementDef)

     j = INDEX( ElementDef(1:n), 'n:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l
       Solver_Def_Dofs(ind,:,1) = l
       IF ( Def_Dofs_Update ) Def_Dofs(ind,1) = MAX(Def_Dofs(ind,1), l)
     END IF

     j = INDEX( ElementDef(1:n), 'e:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l
       Solver_Def_Dofs(ind,:,2) = l
       IF ( Def_Dofs_Update ) Def_Dofs(ind,2) = MAX(Def_Dofs(ind,2), l )
     END IF

     j = INDEX( ElementDef(1:n), 'f:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l
       Solver_Def_Dofs(ind,:,3) = l
       IF ( Def_Dofs_Update ) Def_Dofs(ind,3) = MAX(Def_Dofs(ind,3), l )
     END IF

     j = INDEX( ElementDef(1:n), 'd:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l

       ! Zero value triggers discontinuous approximation within LoadMesh2,
       ! substitute the default negative initialization value to avoid troubles:
       IF (l == 0) l = -1

       Solver_Def_Dofs(ind,:,4) = l
       IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4), l )
     ELSE 
       IF (DG) THEN
         Solver_Def_Dofs(ind,:,4) = 0
         IF ( Def_Dofs_Update ) Def_Dofs(ind,4) = MAX(Def_Dofs(ind,4),0 )
       END IF
     END IF

     j = INDEX( ElementDef(1:n), 'b:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l
       Solver_Def_Dofs(ind,:,5) = l
       IF ( Def_Dofs_Update ) Def_Dofs(ind,5) = MAX(Def_Dofs(ind,5), l )
     END IF

     j = INDEX( ElementDef(1:n), 'p:' )
     IF ( j>0 ) THEN
       IF ( ElementDef(j+2:j+2)=='%' ) THEN
         ! Seeing a p-element definition starting as p:% means that a 
         ! a special keyword construct is used so that the degree of
         ! approximation can be evaluated by calling a MATC function.
         ! This special case is handled elsewhere and we now postpone
         ! setting the right value.
         Solver_Def_Dofs(ind,:,6) = 0
       ELSE
         READ( ElementDef(j+2:), * ) l
         Solver_Def_Dofs(ind,:,6) = l
         IF ( Def_Dofs_Update ) Def_Dofs(ind,6) = MAX(Def_Dofs(ind,6), l )
       END IF
     END IF

!------------------------------------------------------------------------------
   END SUBROUTINE GetDefs
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
! There is no need for calling this unless the element definition is given in
! an equation section or in a body section, or a matc function is used to evaluate
! the order of p-basis, since otherwise the subroutine GetDefs in ModelDescription
! has done the necessary work.
! TO DO: Have just one subroutine for writing def_dofs arrays ?
!------------------------------------------------------------------------------
   SUBROUTINE GetMaxDefs(Model, Mesh, Element, ElementDef, SolverId, BodyId, Def_Dofs)
!------------------------------------------------------------------------------
     TYPE(Model_t) :: Model
     TYPE(MEsh_t) :: Mesh
     TYPE(Element_t) :: Element
     CHARACTER(MAX_NAME_LEN) :: ElementDef
     INTEGER :: SolverId, BodyId, Def_Dofs(:,:)

     INTEGER :: i, j, k, l, n
     INTEGER, POINTER :: ind(:)
     INTEGER, TARGET :: Family(10)
     LOGICAL  :: stat
     REAL(KIND=dp) :: x,y,z
     TYPE(Solver_t), POINTER  :: Solver
     CHARACTER(MAX_NAME_LEN) :: str, ElementDef0


     CALL Info('GetMaxDefs','Checking for other constructs of element definitions', Level=20)

     Family = [1,2,3,4,5,6,7,8,9,10]
     
     Solver => Model % Solvers(SolverId)

     IF ( .NOT. ALLOCATED(Solver % Def_Dofs) ) THEN
       ALLOCATE(Solver % Def_Dofs(10,Model % NumberOfBodies,6))
       Solver % Def_Dofs=-1
       Solver % Def_Dofs(:,:,1)=1
     END IF
     

     ElementDef0 = ElementDef
     DO WHILE(.TRUE.)
       k = INDEX( ElementDef0, '-' )
       IF (k == 1) THEN
         ElementDef0 = ElementDef0(2:)
         k = INDEX( ElementDef0, '-' )
       END IF
         
       IF (k>0) THEN
         !
         ! Read the element definition up to the next flag which specifies the
         ! target element set
         !
         ElementDef = ElementDef0(1:k-1)
       ELSE
         ElementDef = ElementDef0
       END IF

       ! The default assumption is that the given element definition is applied 
       ! to all basic element families (note that the element sets 9 and 10 are
       ! not included since the explicit choice of the target family is 
       ! a part of the element definition string when the target index is
       ! deduced to be 9 or 10).
       !
       ind => Family(1:8)
       !
       ! If the element family is specified, change the target family 
       !       
       IF (SEQL(ElementDef, 'point') )     ind => Family(1:1)
       IF (SEQL(ElementDef, 'line') )      ind => Family(2:2)
       IF (SEQL(ElementDef, 'tri') )       ind => Family(3:3)
       IF (SEQL(ElementDef, 'quad') )      ind => Family(4:4)
       IF (SEQL(ElementDef, 'tetra') )     ind => Family(5:5)
       IF (SEQL(ElementDef, 'pyramid') )   ind => Family(6:6)
       IF (SEQL(ElementDef, 'prism') )     ind => Family(7:7)
       IF (SEQL(ElementDef, 'brick') )     ind => Family(8:8)
       IF (SEQL(ElementDef, 'tri_face') )  ind => Family(9:9)
       IF (SEQL(ElementDef, 'quad_face') ) ind => Family(10:10)

       
       j = INDEX( ElementDef, 'n:' )
       IF ( j>0 ) THEN
         READ( ElementDef(j+2:), * ) l
         Solver % Def_Dofs(ind,BodyId,1) = l
         Def_Dofs(:,1) = MAX(Def_Dofs(:,1), l)
       END IF

       j = INDEX( ElementDef, 'e:' )
       IF ( j>0 ) THEN
         READ( ElementDef(j+2:), * ) l
         Solver % Def_Dofs(ind,BodyId,2) = l
         Def_Dofs(1:8,2) = MAX(Def_Dofs(1:8,2), l )
       END IF

       j = INDEX( ElementDef, 'f:' )
       IF ( j>0 ) THEN
         READ( ElementDef(j+2:), * ) l
         Solver % Def_Dofs(ind,BodyId,3) = l
         Def_Dofs(1:8,3) = MAX(Def_Dofs(1:8,3), l )
       END IF

       j = INDEX( ElementDef, 'd:' )
       IF ( j>0 ) THEN
         READ( ElementDef(j+2:), * ) l

         ! Zero value triggers discontinuous approximation,
         ! substitute the default negative initialization value to avoid troubles:
         IF (l == 0) l = -1

         Solver % Def_Dofs(ind,BodyId,4) = l
         Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4), l )
       ELSE 
         IF ( ListGetLogical( Solver % Values, &
             'Discontinuous Galerkin', stat ) ) THEN
           Solver % Def_Dofs(ind,BodyId,4) = 0
           Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4),0 )
         END IF
       END IF

       j = INDEX( ElementDef, 'b:' )
       IF ( j>0 ) THEN
         READ( ElementDef(j+2:), * ) l
         Solver % Def_Dofs(ind,BodyId,5) = l
         Def_Dofs(1:8,5) = MAX(Def_Dofs(1:8,5), l )
       END IF

       j = INDEX( ElementDef, 'p:' )
       IF ( j>0 ) THEN
         IF ( ElementDef(j+2:j+2) == '%' ) THEN
           n = Element % TYPE % NumberOfNodes
           x = SUM(Mesh % Nodes % x(Element % NodeIndexes))/n
           y = SUM(Mesh % Nodes % y(Element % NodeIndexes))/n
           z = SUM(Mesh % Nodes % z(Element % NodeIndexes))/n
           !          WRITE( str, * ) 'cx= ',i2s(Element % ElementIndex),x,y,z
           str = TRIM(ElementDef(j+3:))//'(cx)'
           x = GetMatcReal(str,4,[1._dp*Element % BodyId,x,y,z],'cx')
           Def_Dofs(1:8,6)  = MAX(Def_Dofs(1:8,6),NINT(x))
           Solver % Def_Dofs(ind,BodyId,6) = NINT(x)
         ELSE
           READ( ElementDef(j+2:), * ) l
           Solver % Def_Dofs(ind,BodyId,6) = l
           Def_Dofs(1:8,6) = MAX(Def_Dofs(1:8,6), l )
         END IF
       END IF

       IF(k>0) THEN
         ElementDef0 = ElementDef0(k+1:)
       ELSE
         EXIT
       END IF
     END DO
!------------------------------------------------------------------------------
  END SUBROUTINE GetMaxDefs
!------------------------------------------------------------------------------




  ! Mark nodes that are associated with at least some boundary element.
  !------------------------------------------------------------------------------
  SUBROUTINE MarkBCNodes(Mesh,BCNode,NoBCNodes)
    TYPE(Mesh_t) :: Mesh
    LOGICAL, ALLOCATABLE :: BCNode(:)
    INTEGER :: NoBCNodes

    INTEGER :: elem
    TYPE(Element_t), POINTER :: Element

    CALL Info('MarkInterfaceNodes','Marking interface nodes',Level=8)

    IF(.NOT. ALLOCATED( BCNode ) ) THEN
      ALLOCATE( BCNode( Mesh % NumberOfNodes ) )
    END IF
    BCNode = .FALSE. 

    DO elem=Mesh % NumberOfBulkElements + 1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      Element => Mesh % Elements( elem )         
      !IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE

      BCNode(Element % NodeIndexes) = .TRUE.
    END DO

    NoBCNodes = COUNT( BCNode )

    CALL Info('MarkBCNodes','Number of BC nodes: '//I2S(NoBCNodes),Level=8)

  END SUBROUTINE MarkBCNodes
!------------------------------------------------------------------------------

  

!------------------------------------------------------------------------------
!> Create a discontinuous mesh over requested boundaries.
!> The nodes are duplicated in order to facilitate the discontinuity.
!> The duplicate nodes are not created by default if the connectivity 
!> of the nodes is needed by other bulk elements than those directly 
!> associated with the discontinuous boundaries. 
!------------------------------------------------------------------------------
 SUBROUTINE CreateDiscontMesh( Model, Mesh, DoAlways )

   TYPE(Model_t) :: Model
   TYPE(Mesh_t) :: Mesh
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
   LOGICAL :: CheckForHalo
   LOGICAL, POINTER :: HaloNode(:)
   TYPE(ValueList_t), POINTER :: BCList
   LOGICAL :: DoneThisAlready = .FALSE.
   CHARACTER(:), ALLOCATABLE :: DiscontFlag
   CHARACTER(*), PARAMETER :: Caller = 'CreateDiscontMesh'

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
   
   CALL Info(Caller,'Creating discontinuous boundaries')

   IF( ActiveBCs > 1 ) THEN
     CALL Warn(Caller,'Be careful when using more than one > Discontinuous Boundary < !')
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
     CALL Warn(Caller,'Missing '//I2S(NoMissingElems)// &
     ' parent elements in partition '//I2S(ParEnv % MyPe)) 
   END IF

   ! Calculate the number of discontinuous nodes and the number of bulk elements 
   ! associated to them. 
   NoDisContElems = COUNT( DiscontElem )
   NoDisContNodes = COUNT( DisContNode ) 
   CALL Info(Caller,'Number of discontinuous boundary elements: '&
       //I2S(NoDisContElems),Level=7)
   CALL Info(Caller,'Number of candicate nodes: '&
       //I2S(NoDisContNodes),Level=7)

   CALL NonGreedyDiscontinuity()
   
   i = ParallelReduction( NoDiscontNodes ) 
   CALL Info(Caller,'Number of discontinuous nodes: '&
       //I2S(i),Level=7)

   IF( i == 0 ) THEN
     CALL Warn(Caller,'Nothing to create, exiting...')
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
   ! boundary. Note that this currently only works in serial!
   IF(.NOT. UseTargetBodies ) THEN
     IF( ParEnv % PEs > 1 ) THEN
       CALL Fatal(Caller,'Please give > Discontinuous Target Bodies < on the BC!')
     END IF
     
     CALL Info(Caller,'Trying to find a dominating parent body',Level=12)

     CandA = -1
     CandB = -1
     DO t=1, NoBoundElems
       IF(.NOT. DisContElem(t) ) CYCLE
       Element => Mesh % Elements(NoBulkElems + t)

       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
         CALL Fatal(Caller,'Alternative strategy requires all parent elements!')
       END IF
       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
         CALL Fatal(Caller,'Alternative strategy requires all parent elements!')
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
     ! This eliminates at the same time the unsuccessful case of zero.
     TargetBody(1) = MAX( CandA, CandB )

     IF( TargetBody(1) > 0 ) THEN
       CALL Info(Caller,&
           'There seems to be a consistent discontinuous body: '&
           //I2S(TargetBody(1)),Level=8)
       UseConsistantBody = .TRUE.
       TargetBodies => TargetBody
     ELSE
       CALL Fatal(Caller,&
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

   CALL Info(Caller,'Number of bulk elements moving: '&
       //I2S(NoMovingElems), Level=8)
   CALL Info(Caller,'Number of bulk elements staying: '&
       //I2S(NoStayingElems), Level=8)

   ! Set discontinuous nodes only if there is a real moving node associated with it
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
         CALL Warn(Caller,'Conflicting moving information')
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
       CALL Fatal(Caller,'Boundary BC has no parents!')
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

   CALL Info(Caller,'Number of related elements moving: '&
       //I2S(NoMovingElems), Level=8 )
   CALL Info(Caller,'Number of related elements staying: '&
       //I2S(NoStayingElems), Level=8 )
   IF( NoUndecided == 0 ) THEN
     CALL Info(Caller,'All elements marked either moving or staying')
   ELSE
     CALL Info(Caller,'Number of related undecided elements: '//I2S(NoUndecided) )
     CALL Warn(Caller,'Could not decide what to do with some boundary elements!')
   END IF


   m = COUNT( DiscontNode .AND. .NOT. MovingNode )
   IF( m > 0 ) THEN
     CALL Info(Caller,'Number of discont nodes not moving: '//I2S(m),Level=12)
   END IF

   m = COUNT( DiscontNode .AND. .NOT. StayingNode )
   IF( m > 0 ) THEN
     CALL Info(Caller,'Number of discont nodes not staying: '//I2S(m),Level=12)
     IF( InfoActive(30) ) THEN
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
   END IF

   !DEALLOCATE( MovingNode, StayingNode )

   ! Now add the new nodes also to the nodes structure
   ! and give the new nodes the same coordinates as the ones
   ! that they were derived from. 
   Mesh % NumberOfNodes = NoNodes + NoDisContNodes   
   CALL EnlargeCoordinates( Mesh ) 

   CALL Info(Caller,'Setting new coordinate positions',Level=12)
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
     CALL Info(Caller,'Creating secondary boundary for Discontinuous gap',Level=10)

     CALL EnlargeBoundaryElements( Mesh, NoDiscontElems ) 

     NoDisContElems = 0
     DO t=1, NoBoundElems

       ! Is this a boundary to be doubled?
       IF(.NOT. DisContElem(t) ) CYCLE

       Element => Mesh % Elements(NoBulkElems + t)
       IF(.NOT. ASSOCIATED(Element) ) THEN
         CALL Fatal(Caller,'Element '//I2S(NoBulkElems+t)//' not associated!')
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
         CALL Fatal(Caller,'Nonzero target boundary must be given for all, if any bc!')
       END IF

       RightElem => Element % BoundaryInfo % Right
       LeftElem => Element % BoundaryInfo % Left 

       NoDisContElems = NoDisContElems + 1              
       j = NoBulkElems + NoBoundElems + NoDisContElems 

       OtherElem => Mesh % Elements( j )
       IF(.NOT. ASSOCIATED(OtherElem) ) THEN
         CALL Fatal(Caller,'Other elem '//I2S(j)//' not associated!')
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

     CALL Info(Caller,'Number of original bulk elements: '&
         //I2S(Mesh % NumberOfBulkElements),Level=10)
     CALL Info(Caller,'Number of original boundary elements: '&
         //I2S(Mesh % NumberOfBoundaryElements),Level=10)
     CALL Info(Caller,'Number of additional boundary elements: '&
         //I2S(NoDisContElems),Level=10)

     Mesh % DiscontMesh = .FALSE.
   ELSE
     Mesh % DisContMesh = .TRUE.
     Mesh % DisContPerm => DisContPerm
     Mesh % DisContNodes = NoDisContNodes 
   END IF

200 CONTINUE

   IF(DoubleBC) THEN
     CALL DropFalseParents()
   END IF
     
   CALL EnlargeParallelInfo(Mesh, DiscontPerm )
   IF( ParEnv % PEs > 1 ) THEN
     m = COUNT( Mesh % ParallelInfo % GlobalDofs == 0) 
     IF( m > 0 ) CALL Warn(Caller,'There are nodes with zero global dof index: '//I2S(m))
   END IF

   IF( DoubleBC .AND. NoDiscontNodes > 0 ) DEALLOCATE( DisContPerm )


   DEALLOCATE( DisContNode, DiscontElem )   
     

 CONTAINS

   ! When indeces change in parents we have to check whether the parents truly are
   ! parents any more!
   !------------------------------------------------------------------------------
   SUBROUTINE DropFalseParents()
     INTEGER :: i,j,t,n,t1,t2,right,hits,nact,npass,nfalse,norphan,torphan
     TYPE(Element_t), POINTER :: Parent, Element
     
     t1 = Mesh % NumberOfBulkElements
     t2 = Mesh % NumberOfBoundaryElements
     nfalse = 0
     norphan = 0
     torphan = 0
     
     DO t = t1+1,t1+t2
       Element => Mesh % Elements(t)
       IF(.NOT. ASSOCIATED(Element % BoundaryInfo) ) CYCLE
       n = Element % TYPE % NumberOfNodes
       nact = 0
       npass = 0
                    
       DO right=0,1
         IF(right==0) THEN
           Parent => Element % BoundaryInfo % Left
         ELSE
           Parent => Element % BoundaryInfo % Right
         END IF
         IF(.NOT. ASSOCIATED(Parent)) CYCLE

         hits = 0
         DO i=1,n
           IF(ANY( Element % NodeIndexes(i) == Parent % NodeIndexes) ) hits = hits + 1
         END DO
         IF( hits == n ) THEN
           nact = nact + 1
           IF(right==1) THEN
             IF(.NOT. ASSOCIATED(Element % BoundaryInfo % Left)) THEN
               Element % BoundaryInfo % Left => Element % BoundaryInfo % Right
               Element % BoundaryInfo % right => NULL()
             END IF
           END IF
         ELSE 
           npass = npass + 1
           IF(right==0) THEN
             Element % BoundaryInfo % Left => NULL()
           ELSE
             Element % BoundaryInfo % Right => NULL()
           END IF
         END IF
       END DO

       IF(npass>0 .AND. nact==0) THEN
         norphan = norphan + 1
         IF( torphan == 0 ) torphan = t
       END IF

       nfalse = nfalse + npass
     END DO

     ! One line for the lot rather than one per element: the count is what
     ! tells whether this is a stray element or the whole boundary.
     IF( norphan > 0 ) THEN
       CALL Warn('DropFalseParents','Boundary elements with no parents of same indexes: '&
           //I2S(norphan)//', first one: '//I2S(torphan))
     END IF

     CALL Info('DropFalseParents','Number of parents no longer parents: '//I2S(nfalse),Level=6)
                      
   END SUBROUTINE DropFalseParents

   
   ! By default all nodes that are associated to elements immediately at the discontinuous 
   ! boundary are treated as discontinuous. However, the user may be not be greedy and release
   ! some nodes from the list that are associated also with other non-discontinuous elements.   
   !-----------------------------------------------------------------------------------------
   SUBROUTINE NonGreedyDiscontinuity()
     INTEGER :: i,i1,i2,j,k
     REAL(KIND=dp) :: Coords(4,3),e1(3),e2(3),phi
     REAL(KIND=dp), ALLOCATABLE :: NodePhi(:)
     INTEGER :: AngleCount(0:36)
     LOGICAL, ALLOCATABLE :: BoundaryNode(:)
     
     IF( NoDiscontNodes == 0 ) RETURN

     ConflictElems = 0

     GreedyBulk = ListGetLogical( Model % Simulation,'Discontinuous Bulk Greedy',Found ) 
     IF(.NOT. Found ) GreedyBulk = .TRUE.     
     
     GreedyBC = ListGetLogical( Model % Simulation,'Discontinuous Boundary Greedy',Found ) 
     IF(.NOT. Found ) GreedyBC = .TRUE.     
          
     IF( .NOT. ( GreedyBC .AND. GreedyBulk ) ) THEN
       CALL Info(Caller,'Applying non-greedy strategies for Discontinuous mesh',Level=12)

       DO t = 1,NoBulkElems+NoBoundElems
         Element => Mesh % Elements(t)
         
         IF( t <= NoBulkElems ) THEN
           IF( GreedyBulk ) CYCLE
           IF( ParentUsed(t) ) CYCLE
         ELSE
           IF( GreedyBC ) CYCLE
           IF( DiscontElem(t-NoBulkElems) ) CYCLE

           ! Check that this is not an external BC
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

       IF( ConflictElems > 0 ) THEN
         CALL Info(Caller,'Conflicting discontinuity in elements: '&
             //I2S(ConflictElems))
       END IF
     END IF

       
     IF( ListGetLogical( Model % Simulation,'Discontinuous Boundary Full Angle',Found ) ) THEN
       CALL Info(Caller,'Computing sum of angles for discontinuous BC',Level=12)

       ALLOCATE(NodePhi(Mesh % NumberOfNodes))
       NodePhi = 0.0_dp

       DO t = 1,NoBoundElems
         Element => Mesh % Elements(NoBulkElems+t)

         IF(.NOT.  DiscontElem(t) ) CYCLE
         
         n = Element % TYPE % ElementCode / 100
         Indexes => Element % NodeIndexes
         Coords(1:n,1) = Mesh % Nodes % y(Indexes(1:n))
         Coords(1:n,2) = Mesh % Nodes % x(Indexes(1:n))
         Coords(1:n,3) = Mesh % Nodes % z(Indexes(1:n))

         DO i = 1, n
           i1 = MODULO(i,n)+1
           i2 = MODULO(n+i-2,n)+1

           e1 = Coords(i1,:)-Coords(i,:)
           e2 = Coords(i2,:)-Coords(i,:)
           
           e1 = e1 / SQRT( SUM( e1**2) )
           e2 = e2 / SQRT( SUM( e2**2) )
           
           ! Cosine angle in radians
           phi = ACOS( SUM( e1 * e2 ) ) 
           
           j = Indexes(i)
           NodePhi(j) = NodePhi(j) + phi
         END DO
       END DO

       ! Move to angles
       NodePhi = 180 * NodePhi / PI
       
       IF( InfoActive(10) ) THEN
         AngleCount = 0
         DO i=1,Mesh % NumberOfNodes
           j = NINT(NodePhi(i)/10)
           AngleCount(j) = AngleCount(j) + 1
         END DO
         DO i=0,36
           j = AngleCount(i)
           IF( j > 0 ) THEN
             CALL Info(Caller,'Angle gat '//I2S(10*i)//' count: '//I2S(j)) 
           END IF
         END DO
       END IF
         
       CALL FindMeshFaces3D(Mesh)
       
       ALLOCATE(BoundaryNode(Mesh % NumberOfNodes) )
       BoundaryNode = .FALSE.
       
       DO t = 1, Mesh % NumberOfFaces
         Element => Mesh % Faces(t)

         i = 0
         IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN           
           IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) i = i+1
           IF( ASSOCIATED( Element % BoundaryInfo % Right) ) i = i+1
         END IF
           
         IF(i==1) THEN
           BoundaryNode(Element % NodeIndexes) = .TRUE.
         END IF
       END DO
       
       i = COUNT( BoundaryNode )
       CALL Info(Caller,'Number of non-internal nodes: '//I2S(i))

       j = 0; k = 0
       DO i = 1, Mesh % NumberOfNodes
         IF(DiscontNode(i) ) THEN
           IF( BoundaryNode(i) ) THEN
             ! On boundary we release the discontinuity when the
             ! angle is ~90 degs i.e. on corner nodes, hopefully. 
             IF( NodePhi(i) < 100.0_dp ) THEN
               DiscontNode(i) = .FALSE.
               j = j+1
             END IF
           ELSE
             ! Elsewhere we release the discontinuity when angle is <360 degs.
             IF( NodePhi(i) < 350.0_dp ) THEN
               DiscontNode(i) = .FALSE.
               k = k+1
             END IF
           END IF
         END IF         
       END DO

       IF(k>0) CALL Info(Caller,'Releasing number of internal boundary nodes: '//I2S(k))
       IF(j>0) CALL Info(Caller,'Releasing number of corner nodes: '//I2S(j))
       
       CALL ReleaseMeshFaceTables( Mesh )
       Mesh % Faces => NULL()

       DEALLOCATE( BoundaryNode, NodePhi ) 

     END IF

     n = NoDiscontNodes
     NoDisContNodes = COUNT( DisContNode ) 

     IF( NoDiscontNodes < n ) THEN
       CALL Info(Caller,'Number of local discontinuous nodes: '&
           //I2S(NoDisContNodes), Level=12)
     ELSE
       CALL Info(Caller,'All candidate nodes used',Level=12)
     END IF

     IF( NoDiscontNodes == 0 ) THEN
       IF( n > 0 .AND. .NOT. GreedyBulk ) THEN
         CALL Info(Caller,'You might want to try the Greedy bulk strategy',Level=3)
       END IF
     END IF

   END SUBROUTINE NonGreedyDiscontinuity  
   
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

   IF(.NOT. ASSOCIATED(Mesh % Nodes % x)) n0 = 0

   pelementsPresent = .FALSE.
   DO i=1,Mesh % NumberOfBulkElements
     IF(isPelement(Mesh % Elements(i))) THEN
       pelementsPresent = .TRUE.; EXIT
     END IF
   END DO

   IF ( Mesh % NumberOfNodes > n0 .OR. n > n0 .AND. pelementsPresent ) THEN
     CALL Info('EnlargeCoordinates','Increasing number of nodes from '&
         //I2S(n0)//' to '//I2S(n),Level=8)

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
       //I2S(n0)//' to '//I2S(n),Level=8)

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
   INTEGER, TARGET :: DiscontPerm(:)

   INTEGER :: nmax,n0,n1,i,j,istat, goffset
   INTEGER, POINTER :: TmpGlobalDofs(:) 
   INTEGER, ALLOCATABLE :: Perm(:)
   LOGICAL, POINTER :: Intf(:)
   TYPE(NeighbourList_t), POINTER :: Nlist(:)

   IF ( ParEnv % PEs <= 1 ) RETURN

   ! As index offset use the number of nodes in the whole mesh
   goffset = ParallelReduction( MAXVAL(Mesh % ParallelInfo % GlobalDofs),2 )

   n0 = SIZE( Mesh % ParallelInfo % GlobalDofs )
   n1 = Mesh % NumberOfNodes 
   IF( n0 >= n1 ) THEN
     CALL Info('EnlargeParallelInfo','No need to grow: '&
         //I2S(n0)//' vs. '//I2S(n1),Level=10)
     RETURN
   END IF
   
   CALL Info('EnlargeParallelInfo','Increasing global numbering size from '&
         //I2S(n0)//' to '//I2S(n1),Level=8)

   ! Create permutation table for the added nodes
   ALLOCATE(Perm(n1)); Perm  = 0
   DO i=1,n0
     IF ( DiscontPerm(i) > 0 ) THEN
       Perm(DiscontPerm(i)+n0) = i
     END IF
   END DO

   ! Create the enlarged set of global nodes indexes
   ALLOCATE( TmpGlobalDofs(n1), STAT=istat )
   IF (istat /= 0) CALL Fatal('EnlargeParallelInfo', 'Unable to allocate TmpGlobalDofs array.')
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
   Intf(1:n0) = Mesh % ParallelInfo % GInterface(1:n0)
   DO i=n0+1,n1
     j = Perm(i)
     IF(j > 0 ) THEN
       Intf(i) = Intf(j) 
     END IF
   END DO
   DEALLOCATE( Mesh % ParallelInfo % GInterface )
   Mesh % ParallelInfo % GInterface => Intf


 END SUBROUTINE EnlargeParallelInfo



 !------------------------------------------------------------------------------
 SUBROUTINE RadiationParallelMeshDistribute(Mesh,nprocs)
 !------------------------------------------------------------------------------
   IMPLICIT NONE

 !------------------------------------------------------------------------------
   TYPE(Mesh_t) :: Mesh
   INTEGER :: nprocs

   INTEGER :: RadiationSurfaces, n_New, n_Coord, n_Coord0, n_Curr, n_NodeInd, max_Coord
   LOGICAL :: Found
   INTEGER :: i,j,k,l,n,ntot,ierr, status(MPI_STATUS_SIZE), narr(ParEnv % PEs)

   REAL(KIND=dp), ALLOCATABLE :: Send_Coord(:), Recv_Coords(:)
   INTEGER, ALLOCATABLE :: ElementNumbers(:), Send_Info(:), Send_Ind(:), &
      Recv_Size(:), Recv_Info(:), Recv_NodeInd(:), Send_Nbr(:), Recv_Nbr(:), cPerm(:)
   LOGICAL, ALLOCATABLE :: CoordsFlag(:)

   TYPE(BoundaryInfo_t), POINTER :: Bi
   TYPE(Element_t), POINTER :: Element, newElements(:)

 !------------------------------------------------------------------------------

   IF(ParEnv % PEs <= 1) RETURN


   ! Get surface elements participating in radiative heat transfer
   ! -------------------------------------------------------------
   ntot = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
   ALLOCATE(ElementNumbers(ntot),CoordsFlag(Mesh % NumberOfNodes))
   ElementNumbers = 0
   CoordsFlag = .FALSE.
   CALL GetMeshRadiationSurfaceInfoA(Mesh,RadiationSurfaces,ElementNumbers,CoordsFlag)

   IF(RadiationSurfaces>0) THEN
     ALLOCATE(Send_Coord(3*COUNT(CoordsFlag)), Send_Nbr(COUNT(CoordsFlag)), &
         cPerm(Mesh % NumberOfNodes+4*ParEnv % PEs*Mesh % NumberOfBoundaryElements))
     cPerm = 0

     ! Extract coordinates of the "owned" radiation elements
     ! -----------------------------------------------------
     n_Coord = 0
     DO i=1,Mesh % NumberOfNodes
       IF ( CoordsFlag(i) ) THEN
         n_Coord = n_Coord + 1
         cPerm(i) = n_Coord
         Send_Coord(3*(n_Coord-1)+1) = Mesh % Nodes % x(i)
         Send_Coord(3*(n_Coord-1)+2) = Mesh % Nodes % y(i)
         Send_Coord(3*(n_Coord-1)+3) = Mesh % Nodes % z(i)
         Send_Nbr(n_Coord) = Mesh % ParallelInfo % GlobalDOFs(i)

           Mesh % ParallelInfo % Ginterface(n_Coord) = .TRUE.

         IF ( ASSOCIATED(Mesh % ParallelInfO % NeighbourList(i) % Neighbours)) THEN
           DEALLOCATE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
         END IF

         ALLOCATE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours(ParEnv % PEs))
         Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % myPE
         k = 1
         DO j=0,ParEnv % PEs-1
           IF ( j==ParEnv % Mype) CYCLE
           k = k + 1 
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours(k) = j
         END DO
       END IF
     END DO

     ! Extract topolgy of the "owned" radiation elements
     ! --------------------------------------------------
     ALLOCATE(Send_Info(3*RadiationSurfaces), Send_ind(4*RadiationSurfaces))
     n_NodeInd  = 0
     DO i=1,RadiationSurfaces
       j = ElementNumbers(i)
       Element  => Mesh % Elements(j)
       n = ElemenT % Type % NumberOfNodes
       Element % PartIndex = ParEnv % myPE

       Send_Info(3*(i-1)+1) = Element % Type % ElementCode
       Send_Info(3*(i-1)+2) = Element % BoundaryInfo % Constraint
       Send_Info(3*(i-1)+3) = Element % GElementIndex
       Send_Ind(n_NodeInd+1:n_NodeInd+n) = cPerm(Element % NodeIndexes)
       n_NodeInd = n_NodeInd + n
     END DO
     CALL CheckBuffer(ParEnv % PEs*(3*RadiationSurfaces+n_NodeInd+8*n_Coord+MPI_BSEND_OVERHEAD))
   ELSE
     CALL CheckBuffer(1024+MPI_BSEND_OVERHEAD) ! just something
   END IF

   ! Distribute the extracted information
   ! ------------------------------------
   DO i=0,nprocs-1
     IF (i==ParEnv % myPE) CYCLE

     CALL MPI_BSEND( RadiationSurfaces,1,MPI_INTEGER,i,12000,ELMER_COMM_WORLD,ierr )
     IF ( RadiationSurfaces>0 ) THEN
       CALL MPI_BSEND( Send_Info,3*RadiationSurfaces,MPI_INTEGER,i,12001,ELMER_COMM_WORLD,ierr )
       CALL MPI_BSEND( Send_Ind, n_NodeInd,MPI_INTEGER,i,12002,ELMER_COMM_WORLD,ierr )
       CALL MPI_BSEND( Send_Nbr,n_Coord,MPI_INTEGER,i,12003,ELMER_COMM_WORLD,ierr )
       CALL MPI_BSEND( Send_Coord,3*n_Coord,MPI_DOUBLE_PRECISION,i,12004,ELMER_COMM_WORLD,ierr )
     END IF
   END DO

   ! Receive element counts from around (if we don't own any, ignore others too...)
   ! ------------------------------------------------------------------------------
   ALLOCATE(Recv_Size(0:ParEnv % PEs-1))
   Recv_Size = 0
   DO i=0,nprocs-1
     IF (i==ParEnv % myPE) CYCLE
     CALL MPI_RECV( Recv_Size(i),1,MPI_INTEGER,i,12000,ELMER_COMM_WORLD,status,ierr )
   END DO

   ! If we are not participating in radiation or none else is, skip the rest
   !-------------------------------------------------------------------------
   IF (RadiationSurfaces==0 .OR. SUM(Recv_Size)==0 ) RETURN

   ! Receive surface elements
   ! ------------------------
   n_Curr  = ntot
   n_New   = SUM(Recv_Size) + n_Curr
   n_Coord = Mesh % NumberOfNodes

   ! Re-allocate mesh structures to contain the received surface elements 
   ! --------------------------------------------------------------------
   BLOCK
     TYPE(NeighbourList_t), POINTER :: x(:)
     LOGICAL, POINTER :: y(:)
     INTEGER, POINTER :: z(:)
     REAL(KIND=dp), POINTER :: xc(:),yc(:),zc(:)

     ALLOCATE( x(n_Coord + 4*SUM(Recv_Size)) )
     x(1:n_Coord) = Mesh % ParallelInfo % NeighbourList
     DO i=n_Coord+1, n_Coord + SUM(Recv_Size)
       x(i) % Neighbours=> Null()
     END DO
     DEALLOCATE(Mesh % ParallelInfo % NeighbourList)
     Mesh % ParallelInfo % NeighbourList => x

     ALLOCATE( y(n_Coord + 4*SUM(Recv_Size)) )
     y(1:n_Coord) = Mesh % ParallelInfo % Ginterface
     DEALLOCATE( Mesh % ParallelInfo % Ginterface )
     Mesh % ParallelInfo % Ginterface => y

     ALLOCATE( z(n_Coord + 4*SUM(Recv_Size)) )
     z(1:n_Coord) = Mesh % ParallelInfo % GlobalDofs
     DEALLOCATE( Mesh % ParallelInfo % GlobalDofs )
     Mesh % ParallelInfo % GlobalDofs => z

     ALLOCATE( xc(n_Coord + SUM(Recv_Size)), &
               yc(n_Coord + SUM(Recv_Size)), &
               zc(n_Coord + SUM(Recv_Size)) )

     xc(1:n_Coord) = Mesh % Nodes % x
     yc(1:n_Coord) = Mesh % Nodes % y
     zc(1:n_Coord) = Mesh % Nodes % z
     DEALLOCATE( Mesh % Nodes % x, Mesh % Nodes % y, Mesh % Nodes % z)
     Mesh % Nodes % x => xc
     Mesh % Nodes % y => yc
     Mesh % Nodes % z => zc
   END BLOCK

   ALLOCATE( newElements(n_New) )
   newElements(1:n_Curr) = Mesh % Elements
   DO i=1,Mesh % NumberOfBoundaryElements
     Element => NewElements(i+Mesh % NumberOfBulkElements)
     Bi => Element % BoundaryInfo
     IF ( ASSOCIATED(Bi) ) THEN
       IF(ASSOCIATED(Bi % Left))  BI % Left  => NewElements(Bi % Left % ElementIndex)
       IF(ASSOCIATED(Bi % Right)) BI % Right => NewElements(Bi % Right % ElementIndex)
     END IF
   END DO
   DEALLOCATE(Mesh % Elements)
   Mesh % Elements => newElements

   ! Receive the elements from other partitions
   ! ------------------------------------------
   n = MAXVAL(Recv_Size)
   ALLOCATE(Recv_Info(3*n), Recv_NodeInd(4*n), Recv_Coords(12*n), Recv_Nbr(4*n))

   n_Coord = Mesh % NumberOfNodes
   DO i=0,nprocs-1
     IF (Recv_Size(i) <= 0) CYCLE

     CALL MPI_RECV( Recv_Info,3*Recv_Size(i),MPI_INTEGER,i,12001,ELMER_COMM_WORLD,status,ierr )
     CALL MPI_RECV( Recv_NodeInd,4*Recv_Size(i),MPI_INTEGER,i,12002,ELMER_COMM_WORLD,status,ierr )

     CALL MPI_RECV( Recv_Nbr,4*Recv_Size(i),MPI_INTEGER,i,12003,ELMER_COMM_WORLD,status,ierr )
     CALL MPI_GET_COUNT(status, MPI_INTEGER, n, ierr )

     CALL MPI_RECV( Recv_Coords,12*Recv_Size(i),MPI_DOUBLE_PRECISION,&
          i,12004,ELMER_COMM_WORLD,status,ierr )

     ! Insert the received elements to the (already mostly re-allocated) mesh strucres
     ! -------------------------------------------------------------------------------
     n_Coord0 = n_Coord

     BLOCK
       INTEGER, ALLOCATABLE :: Gdofs(:), Gorder(:)

       Gdofs = Mesh % ParallelInfo % GlobalDOFs
       Gorder = [(j, j=1,n_Coord0)]
       CALL Sorti(n_Coord0,Gdofs,Gorder)

       ! Insert nodes
       ! ------------
       DO j=1,n

         k = SearchNode(Mesh % ParallelInfo,Recv_Nbr(j),1,n_Coord0,Gorder)
         IF(k>0) THEN
           cPerm(j)=k;

           IF(.NOT.ASSOCIATED(Mesh % ParallelInfo % NeighbourList(k) % Neighbours)) STOP 'a'

           l = SIZE(Mesh % ParallelInfo % NeighbourList(k) % Neighbours)
           narr(1:l) = Mesh % ParallelInfo % NeighbourList(k) % Neighbours

           IF (ALL(narr(1:l) /= i)) THEN
             l = l +1
             narr(l) = i
             DEALLOCATE(Mesh % ParallelInfo % NeighbourList(k) % Neighbours)
             ALLOCATE(Mesh % ParallelInfo % NeighbourList(k) % Neighbours(l))
             Mesh % ParallelInfo % NeighbourList(k) % Neighbours = narr(1:l)
           END IF

           IF (ALL(narr(1:l) /= ParEnv % myPE)) THEN
             l = l +1
             narr(l) = ParEnv % myPE
             DEALLOCATE(Mesh % ParallelInfo % NeighbourList(k) % Neighbours)
             ALLOCATE(Mesh % ParallelInfo % NeighbourList(k) % Neighbours(l))
             Mesh % ParallelInfo % NeighbourList(k) % Neighbours = narr(1:l)
           END IF

           IF ( .NOT. Mesh % ParallelInfo % Ginterface(k)) THEN
             Mesh % ParallelInfo % Ginterface(k) = .TRUE.
             Mesh % ParallelInfo % NumberOfIfDOFs = Mesh % ParallelInfo % NumberOfIfDOFs+1
           END IF
           CYCLE
         END IF

         n_Coord = n_Coord + 1

         cPerm(j)=n_Coord
         Mesh % Nodes % x(n_Coord) = Recv_Coords(3*(j-1)+1)
         Mesh % Nodes % y(n_Coord) = Recv_Coords(3*(j-1)+2)
         Mesh % Nodes % z(n_Coord) = Recv_Coords(3*(j-1)+3)

         Mesh % ParallelInfo % NumberOfIfDOFs = Mesh % ParallelInfo % NumberOfIfDOFs+1
         Mesh % ParallelInfo % Ginterface(n_Coord) = .TRUE.
         Mesh % ParallelInfo % GlobalDofs(n_Coord) = Recv_Nbr(j)

         IF(ASSOCIATED(Mesh % ParallelInfo % NeighbourList(n_Coord) % Neighbours)) STOP 'b'

         ALLOCATE(Mesh % ParallelInfo % NeighbourList(n_Coord) % Neighbours(2))
         Mesh % ParallelInfo % NeighbourList(n_Coord) % Neighbours(1) = i
         Mesh % ParallelInfo % NeighbourList(n_Coord) % Neighbours(2) = ParEnv % myPE
       END DO
     END BLOCK

     ! Insert elements
     ! ---------------
     k = 0
     DO j=1,Recv_Size(i)
       Element => Mesh % Elements(j+n_Curr)

       Element % Type => GetElementType(Recv_Info(3*(j-1)+1))
       n = Element % Type % NumberOfNodes

       ALLOCATE(Element % BoundaryInfo)
       Element % BoundaryInfo % Constraint = Recv_Info(3*(j-1)+2)
       Element % BoundaryInfo % Left => Null()
       Element % BoundaryInfo % Right => Null()

       Element % PartIndex = i
       Element % BodyId = 0
       Element % ElementIndex  = j+n_Curr
       Element % GElementIndex = Recv_Info(3*(j-1)+3)

       ALLOCATE(Element % NodeIndexes(n))
       Element % NodeIndexes = cPerm(Recv_NodeInd(k+1:k+n))
       k = k + n
     END DO
     n_Curr = n_Curr + Recv_Size(i)
   END DO

   Mesh % NumberOFnodes = n_Coord

   ! Try reset the owner of a node (first entry in the node's Neighbours-array)
   ! to some commonly knowable task
   ! ---------------------------------------------------------------------------
   BLOCK
     INTEGER, POINTER :: Neighbours(:)
     INTEGER, ALLOCATABLE :: gCount(:)

     ALLOCATE( gCount(Mesh % NumberOfNodes) ); gCount = 0

     DO i=1,Mesh % NumberOfBoundaryElements+SUM(Recv_Size)
       Element => Mesh % Elements(i+Mesh % NumberOfBulkElements)
       IF ( .NOT. RadiationCheck(Element)) CYCLE
       DO j=1,Element % Type % NumberOfNodes
         n = Element % NodeIndexes(j)
         gCount(n) = MAX(Element % GElementIndex,gCount(n))
       END DO
     END DO

     DO i=1,Mesh % NumberOfBoundaryElements+SUM(Recv_Size)
       Element => Mesh % Elements(i+Mesh % NumberOfBulkElements)
       IF ( .NOT. RadiationCheck(Element)) CYCLE
       DO j=1,Element % Type % NumberOfNodes
         n = Element % NodeIndexes(j)

         IF ( Element % GElementIndex /= gCount(n)) CYCLE

         Neighbours => Mesh % ParallelInfo % NeighbourList(n) % Neighbours
         DO k=1,SIZE(Neighbours)
           IF(Neighbours(k) == Element % PartIndex) EXIT
         END DO
         if ( k>SIZE(Neighbours) ) stop 'fail0'

         l = Neighbours(1); Neighbours(1) = Element % PartIndex; Neighbours(k) = l
         if ( Element % PartIndex == parenv % mype) then
            IF ( .NOT.ASSOCIATED(element % boundaryinfo % left ) ) stop 'fail1'
            IF ( neighbours(1) /= parenv % mype ) stop 'fail2'
         end if
       END DO
     END DO
   END BLOCK

   Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements + SUM(Recv_Size)

CONTAINS

 !------------------------------------------------------------------------------
 SUBROUTINE GetMeshRadiationSurfaceInfoA(Mesh,RadiationSurfaces,ElementNumbers,CoordsFlag)
 !------------------------------------------------------------------------------
   IMPLICIT NONE

   TYPE(ValueList_t), POINTER :: BC
   INTEGER ::  ElementNumbers(:)
   LOGICAL ::  CoordsFlag(:)
   TYPE(Mesh_t) :: Mesh
   INTEGER :: i,j,t,n, RadiationSurfaces, nbulk
   LOGICAL :: Found
   TYPE(Element_t), POINTER :: Element

   nBulk = Mesh % NumberOfBulkElements
   RadiationSurfaces = 0
   ElementNumbers = 0
   CoordsFlag = .FALSE.
   DO i=1,Mesh % NumberOfBoundaryElements
     Element => Mesh % Elements(nBulk+i)
     IF (RadiationCheck(Element)) THEN
       RadiationSurfaces = RadiationSurfaces + 1
       CoordsFlag(Element % NodeIndexes) = .TRUE.
       ElementNumbers(RadiationSurfaces) = i + nBulk
     END IF
   END DO
 !------------------------------------------------------------------------------
 END SUBROUTINE GetMeshRadiationSurfaceInfoA
 !------------------------------------------------------------------------------

 !------------------------------------------------------------------------------
 END SUBROUTINE RadiationParallelMeshDistribute
 !------------------------------------------------------------------------------

 
 !------------------------------------------------------------------------------
 !> Function to load mesh from disk.
 !------------------------------------------------------------------------------
 !------------------------------------------------------------------------------


 !> Prepare a clean nodal mesh as it comes after being loaded from disk.
 !> Study the non-nodal elements (face, edge, DG, and p-elements)
 !> Create parallel info for the non-nodal elements
 !> Enlarge the coordinate vectors for p-elements.
 !> Generate static projector for periodic BCS.
 !-------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Transfer coordinate and time from one mesh toanother when swapping meshes
!> for some reason.
!------------------------------------------------------------------------------  
  SUBROUTINE TransferCoordAndTime(M1,M2)
    TYPE(Solver_t), POINTER :: Solver => Null()
    TYPE(Mesh_t) :: M1,M2
    TYPE(Variable_t), POINTER :: DtVar, V

     CALL VariableAdd( M2 % Variables, M2,Solver, &
           'Coordinate 1',1, M2 % Nodes % x )

     CALL VariableAdd(M2 % Variables,M2,Solver, &
           'Coordinate 2',1, M2 % Nodes % y )

     CALL VariableAdd(M2 % Variables,M2,Solver, &
          'Coordinate 3',1,M2 % Nodes % z )

     V => VariableGet( M1 % Variables, 'Time' )     
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Time', 1, V % Values )

     V => VariableGet( M1 % Variables, 'Periodic Time' )
     IF( ASSOCIATED( V ) ) THEN
       CALL VariableAdd( M2 % Variables, M2, Solver, 'Periodic Time', 1, V % Values)
     END IF
     V => VariableGet( M1 % Variables, 'Periodic Cycle' )
     IF( ASSOCIATED( V ) ) THEN
       CALL VariableAdd( M2 % Variables, M2, Solver, 'Periodic Cycle', 1, V % Values)
     END IF
       
     V => VariableGet( M1 % Variables, 'Timestep' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Timestep', 1, V % Values )

     V => VariableGet( M1 % Variables, 'Timestep size' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Timestep size', 1, V % Values )

     V => VariableGet( M1 % Variables, 'Timestep interval' )
     CALL VariableAdd( M2 % Variables, M2, Solver, 'Timestep interval', 1, V % Values )

     ! Save some previous timesteps for variable timestep multistep methods
     V => VariableGet( M1 % Variables, 'Timestep size' )
     DtVar => VariableGet( M2 % Variables, 'Timestep size' )
     DtVar % PrevValues => V % PrevValues

     V => VariableGet( M1 % Variables, 'nonlin iter' )
     CALL VariableAdd( M2 % Variables, M2, Solver, &
         'nonlin iter', 1, V % Values )
     
     V => VariableGet( M1 % Variables, 'coupled iter' )
     CALL VariableAdd( M2 % Variables, M2, Solver, &
         'coupled iter', 1, V % Values )
     
     V => VariableGet( M1 % Variables, 'partition' )
     IF( ASSOCIATED( V ) ) THEN
       CALL VariableAdd( M2 % Variables, M2, Solver, 'Partition', 1, V % Values )
     END IF
     
     V => VariableGet( M1 % Variables, 'scan' )
     IF( ASSOCIATED( V ) ) THEN
       CALL VariableAdd( M2 % Variables, M2, Solver, 'scan', 1, V % Values)
     END IF
     V => VariableGet( M1 % Variables, 'finish' )
     IF( ASSOCIATED( V ) ) THEN
       CALL VariableAdd( M2 % Variables, M2, Solver, 'finish', 1, V % Values)
     END IF
     V => VariableGet( M1 % Variables, 'produce' )
     IF( ASSOCIATED( V ) ) THEN
       CALL VariableAdd( M2 % Variables, M2, Solver, 'produce', 1, V % Values)
     END IF
     V => VariableGet( M1 % Variables, 'run' )
     IF( ASSOCIATED( V ) ) THEN
       CALL VariableAdd( M2 % Variables, M2, Solver, 'run', 1, V % Values)
     END IF
     
!------------------------------------------------------------------------------
   END SUBROUTINE TransferCoordAndTime
!------------------------------------------------------------------------------


  !-------------------------------------------------------------------------------
  !> Communicate logical tag related to mesh or linear system.
  !> This could related to setting Neumann BCs to zero, for example.
  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateParallelSystemTag(ParallelInfo,Ltag,Itag,ParOper)
  !-------------------------------------------------------------------------------
     TYPE (ParallelInfo_t), POINTER :: ParallelInfo
     LOGICAL, OPTIONAL :: LTag(:)   !< Logical tag, if used
     INTEGER, OPTIONAL :: ITag(:)   !< Integer tag, if used
     INTEGER, OPTIONAL :: ParOper            !< If integer tag is used, we can also have an operator

     LOGICAL, POINTER :: IsNeighbour(:)
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:), s_i(:,:), r_i(:)
     INTEGER :: i,j,k,l,n,nn,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)
     INTEGER :: NewZeros, nsize
     LOGICAL :: UseL, GotIt
     INTEGER :: CommI
     
     IF( ParEnv % PEs<=1 ) RETURN
   
     UseL = PRESENT(LTag)
     IF(.NOT. (UseL .NEQV. PRESENT(Itag)) ) THEN
       CALL Fatal('CommunicateParallelSystemTag','Give either logical or integer tag!')
     END IF
     CommI = -1
     IF(.NOT. UseL) THEN
       IF(PRESENT(ParOper)) CommI = ParOper
     END IF
     
     ! The caller may hand us the ParallelInfo of a matrix that was created
     ! but never passed through ParallelInitMatrix. Both the structure itself
     ! and GInterface are pointers, so SIZE() below would segfault rather than
     ! say anything useful about which matrix is at fault.
     IF(.NOT. ASSOCIATED( ParallelInfo ) ) THEN
       CALL Fatal('CommunicateParallelSystemTag',&
           'ParallelInfo not created, call ParallelInitMatrix for the matrix first!')
     END IF
     IF(.NOT. ASSOCIATED( ParallelInfo % GInterface ) ) THEN
       CALL Fatal('CommunicateParallelSystemTag',&
           'ParallelInfo % GInterface not created, call ParallelInitMatrix for the matrix first!')
     END IF

     nsize = SIZE( ParallelInfo % GInterface)
     IF( PRESENT(Ltag) ) THEN
       nsize = MIN(nsize, SIZE(Ltag) )
     ELSE
       nsize = MIN(nsize, SIZE(Itag) )
     END IF
     
     ALLOCATE( fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )
     
     ! Mark the neighbouring entities
     IF(ASSOCIATED( ParEnv % IsNeighbour ) ) THEN
       IsNeighbour => ParEnv % IsNeighbour
     ELSE
       ! We may want to call this even though neighbours have not been set
       ALLOCATE( IsNeighbour(ParEnv % PEs) )
       IsNeighbour = .FALSE.
       DO i=1,nsize 
         DO j=1,SIZE(ParallelInfo % Neighbourlist(i) % Neighbours)
           k = ParallelInfo % Neighbourlist(i) % Neighbours(j)
           IF ( k == ParEnv % MyPE ) CYCLE
           IsNeighbour(k+1) = .TRUE.
         END DO
       END DO
     END IF
     
     nn = 0
     ineigh = 0
     DO i=0, ParEnv % PEs-1
       k = i+1
       IF(.NOT.ParEnv % Active(k) ) CYCLE
       IF(i == ParEnv % myPE) CYCLE
       IF(.NOT. IsNeighbour(k) ) CYCLE
       nn = nn + 1
       fneigh(nn) = k
       ineigh(k) = nn
     END DO

     IF(.NOT. ASSOCIATED( ParEnv % IsNeighbour ) ) THEN
       DEALLOCATE(IsNeighbour)
     END IF
     
     ! Count the maximum number of enties to sent 
     IF( UseL ) THEN
       n = COUNT(LTag(1:nsize) .AND. ParallelInfo % GInterface(1:nsize))
     ELSE
       n = COUNT((ITag(1:nsize) /= 0) .AND. ParallelInfo % GInterface(1:nsize))
     END IF

     ! Allocate for the data to sent (s_e) and receive (r_e)
     ALLOCATE( s_e(n, nn ), r_e(n) )
     s_e = 0
     IF( CommI >= 0 ) THEN
       ALLOCATE( s_i(n, nn), r_i(n) )
       s_i = 0
     END IF

     IF( CommI >= 0) THEN
       CALL CheckBuffer( nn*6*n )
     ELSE
       CALL CheckBuffer( nn*3*n )
     END IF
       
     ii = 0
     DO i=1, nsize
       IF( UseL ) THEN
         GotIt = LTag(i) .AND. ParallelInfo % GInterface(i)
       ELSE
         GotIt = Itag(i) /= 0 .AND. ParallelInfo % GInterface(i)
       END IF
       IF(.NOT. GotIt) CYCLE
       
       DO j=1,SIZE(ParallelInfo % Neighbourlist(i) % Neighbours)
         k = ParallelInfo % Neighbourlist(i) % Neighbours(j)
         IF ( k == ParEnv % MyPE ) CYCLE
         k = k + 1
         k = ineigh(k)
         IF ( k> 0) THEN
           ii(k) = ii(k) + 1
           s_e(ii(k),k) = ParallelInfo % GlobalDOFs(i)
           IF( CommI >= 0 ) THEN
             s_i(ii(k),k) = Itag(i)
           END IF
         END IF
       END DO
     END DO

     DO i=1, nn
       j = fneigh(i) 
       ! Sent size data
       CALL MPI_BSEND( ii(i),1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD,ierr )
       IF( ii(i) > 0 ) THEN
         ! Sent the global index 
         CALL MPI_BSEND( s_e(1:ii(i),i),ii(i),MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,ierr )
         IF( CommI >= 0 ) THEN
           ! Sent the value of the integer tag, if requested
           CALL MPI_BSEND( s_i(1:ii(i),i),ii(i),MPI_INTEGER,j-1,112,ELMER_COMM_WORLD,ierr )
         END IF
       END IF
     END DO

     NewZeros = 0
     
     DO i=1, nn
       j = fneigh(i)
       ! Receive size of data coming from partition "j"
       CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
       IF ( n>0 ) THEN
         IF( n>SIZE(r_e)) THEN
           DEALLOCATE(r_e)
           ALLOCATE(r_e(n))
           IF(CommI >= 0) THEN
             DEALLOCATE(r_i)
             ALLOCATE(r_i(n))
           END IF
         END IF

         ! Receive the global index
         CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
         IF( CommI >= 0) THEN
           ! Receive the value of the integer tag, if requested
           CALL MPI_RECV( r_i,n,MPI_INTEGER,j-1,112,ELMER_COMM_WORLD,status,ierr )
         END IF
         DO j=1,n
           ! Check that the entry exists in the matrix
           k = SearchNode( ParallelInfo, r_e(j), Order=ParallelInfo % Gorder )
           IF ( k>0 ) THEN
             IF( UseL ) THEN
               IF(.NOT. LTag(k)) THEN
                 LTag(k) = .TRUE.
                 NewZeros = NewZeros + 1
               END IF
             ELSE
               IF( CommI == 0 ) THEN
                 ! "i" here is the neighbour loop counter, the row to update
                 ! is the one SearchNode found, i.e. "k".
                 Itag(k) = Itag(k) + r_i(j)
               ELSE IF( CommI == 1 ) THEN
                 ITag(k) = MIN(r_i(j),Itag(k))
               ELSE IF( CommI == 2 ) THEN
                 ITag(k) = MAX(r_i(j),Itag(k))
               ELSE IF( Itag(k) == 0 ) THEN
                 ITag(k) = 1
               END IF
               NewZeros = NewZeros + 1
             END IF
           END IF
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e )
     IF(CommI >= 0) DEALLOCATE(s_i, r_i)

     !PRINT *,'New Zeros:',ParEnv % MyPe, NewZeros
     
  !-------------------------------------------------------------------------------
   END SUBROUTINE CommunicateParallelSystemTag
  !-------------------------------------------------------------------------------

 

 ! This subroutine fixes the global indexing of the mesh when the same mesh has been loaded to the
 ! for multiple partitions.
 !-------------------------------------------------------------------------------------------------
 SUBROUTINE SetMeshPartitionOffset(Mesh,nParMesh)
   TYPE(Mesh_t), TARGET :: Mesh
   INTEGER :: nParMesh
   
   INTEGER :: Offset
   INTEGER :: i,n,ierr,iParExt,nParExt
   TYPE(ParallelInfo_t), POINTER :: PI

   CALL Info('SetMeshPartitionOffset','Setting offset when same mesh loaded for multiple partitions!')
   
   IF( nParMesh < 1 .OR. nParMesh >= ParEnv % PEs ) THEN
     CALL Fatal('SetMeshPartitionOffset','Invalid value of parameter nParMesh: '//I2S(nParMesh))
   END IF
   IF( MODULO(ParEnv % PEs, nParMesh ) /= 0 ) THEN
     CALL Fatal('SetMeshPartitionOffset','Number of partitions should be divisible with: '//I2S(nParMesh))
   END IF
   
   nParExt = ParEnv % PEs / nParMesh
   iParExt = ParEnv % MyPe / nParMesh

   
   PI => Mesh % ParallelInfo
   
   ! update neighbourist for partitions with an offset   
   DO i=1,Mesh % NumberOfNodes 
     IF (ASSOCIATED(PI % NeighbourList(i) % Neighbours)) THEN
       PI % NeighbourList(i) % Neighbours = &
           PI % NeighbourList(i) % Neighbours + iParExt * nParMesh
     END IF
   END DO
 
   ! Set offset for global node indexes, first find the max node index and then add the offset
   i = MAXVAL(PI % GlobalDofs )                
   CALL MPI_ALLREDUCE(i,n,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
   DO i=1,Mesh % NumberOfNodes
     PI % GlobalDofs(i) = PI % GlobalDofs(i) + iParExt * n
   END DO
   
   ! Set offset for global element indexes, first find the max element index the add the offset   
   i = MAXVAL(Mesh % Elements(:) % GElementIndex )  
   CALL MPI_ALLREDUCE(i,n,1, MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)   
   DO i=1,Mesh % NumberOfBulkElements
     Mesh % Elements(i) % GElementIndex = Mesh % Elements(i) % GElementIndex + iParExt * n
     Mesh % Elements(i) % PartIndex = Mesh % Elements(i) % PartIndex + iParExt * nParMesh
   END DO
   
 END SUBROUTINE SetMeshPartitionOffset
   
 
!------------------------------------------------------------------------------
  SUBROUTINE SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs,NeedEdges)
!------------------------------------------------------------------------------
    INTEGER, OPTIONAL :: EdgeDOFs(:), FaceDOFs(:)
    TYPE(Mesh_t) :: Mesh
    INTEGER, OPTIONAL :: indofs(:,:)
    LOGICAL, OPTIONAL :: NeedEdges
!------------------------------------------------------------------------------
    INTEGER :: i,j,el_id
    TYPE(Element_t), POINTER :: Element, Edge, Face
    LOGICAL :: AssignEdges, pAlloc
!------------------------------------------------------------------------------

    CALL FindMeshEdges(Mesh)
    
    AssignEdges = .FALSE.
    IF (PRESENT(NeedEdges)) AssignEdges = NeedEdges
    
    CALL Info('SetMeshEdgeFaceDofs','Setting edge and face dofs for elements!',Level=20)
    
    DO i=1,Mesh % NumberOFBulkElements
       Element => Mesh % Elements(i)
       
       IF(ASSOCIATED(Element % EdgeIndexes)) THEN
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
               ELSE IF(ASSOCIATED(Edge % BoundaryInfo % Right)) THEN
                 CALL AssignLocalNumber(Edge, Edge % BoundaryInfo % Right, Mesh)
               END IF
             END IF
            ! Other element types, which need edge dofs
            IF(PRESENT(EdgeDOFs)) THEN
              Edge % BDOFs = MAX(EdgeDOFs(i), Edge % BDOFs)
            ELSE
              Edge % BDOFs = Max(1, Edge % BDOFs)
            END IF

            ! Get maximum dof for edges
            Mesh % MinEdgeDOFs = MIN(Edge % BDOFs, Mesh % MinEdgeDOFs)
            Mesh % MaxEdgeDOFs = MAX(Edge % BDOFs, Mesh % MaxEdgeDOFs)
         END DO
       END IF
       IF ( Mesh % MinEdgeDOFs > Mesh % MaxEdgeDOFs ) Mesh % MinEdgeDOFs = MEsh % MaxEdgeDOFs

       ! Iterate each face of element
       IF(.NOT. ASSOCIATED(Element % FaceIndexes)) CYCLE

       DO j=1,Element % TYPE % NumberOfFaces
          Face => Mesh % Faces( Element % FaceIndexes(j) )

          IF(ANY(Face % EdgeIndexes==0)) CYCLE

          ! Set attributes of p element faces
          IF ( ASSOCIATED(Element % PDefs) ) THEN
             ! Set face polynomial degree and dofs
             Face % PDefs % P = MAX(Element % PDefs % P, Face % PDefs % P)
             ! Get number of face dofs
             Face % BDOFs = MAX(Face % BDOFs, getFaceDOFs(Element, Face % PDefs % P, j,Face) )
             Face % PDefs % isEdge = .TRUE.
             Face % PDefs % GaussPoints = getNumberOfGaussPointsFace( Face, Mesh )
             IF (ASSOCIATED(Face % BoundaryInfo % Left) ) THEN
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Left, Mesh)
             ELSE
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Right, Mesh)
             END IF
           END IF
           IF (PRESENT(FaceDOFs)) THEN
             !
             ! NOTE: This depends on what dofs have been introduced
             ! by using the construct "-quad_face b: ..." and
             ! "-tri_face b: ..."
             !
             el_id = face % TYPE % ElementCode / 100
             Face % BDOFs = MAX(FaceDOFs(i), Face % BDOFs)
             IF ( PRESENT(inDOFs) ) Face % BDOFs = MAX(Face % BDOFs, InDOFs(el_id+6,5))
          END IF
             
          ! Get maximum dof for faces
          Mesh % MinFaceDOFs = MIN(Face % BDOFs, Mesh % MinFaceDOFs)
          Mesh % MaxFaceDOFs = MAX(Face % BDOFs, Mesh % MaxFaceDOFs)
       END DO
    END DO
    IF ( Mesh % MinFaceDOFs > Mesh % MaxFaceDOFs ) Mesh % MinFaceDOFs = Mesh % MaxFaceDOFs

    ! Set local edges for boundary elements

    CALL Info('SetMeshEdgeFaceDofs','Setting local edges for boundary elements',Level=20)

    DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
      
       ! Here set local number and copy attributes to this boundary element for left parent.
       pAlloc = .FALSE.
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
         ! Local edges are only assigned for p elements
         IF (ASSOCIATED(Element % BoundaryInfo % Left % PDefs)) THEN
           pAlloc = .TRUE.
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
           IF(.NOT. pAlloc) THEN
             CALL AllocatePDefinitions(Element)
             Element % PDefs % isEdge = .TRUE.
             CALL AssignLocalNumber(Element, Element % BoundaryInfo % Right, Mesh)
           END IF
         END IF
       END IF

       IF (AssignEdges) THEN
         IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
           CALL AssignLocalNumber(Element,Element % BoundaryInfo % Left, Mesh, NoPE=.TRUE.)
         END IF
         IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
           CALL AssignLocalNumber(Element,Element % BoundaryInfo % Right, Mesh, NoPE=.TRUE.)
         END IF
       END IF
     END DO

    CALL Info('SetMeshEdgeFaceDofs','All done',Level=25)

     
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

   DO i=1,Mesh % NumberOfBulkElements
     Element => Mesh % Elements(i)

     ! Set gauss points for each p element
     IF ( ASSOCIATED(Element % PDefs) ) THEN
       Element % PDefs % GaussPoints = getNumberOfGaussPoints( Element, Mesh )
     END IF

     Mesh % MaxBDOFs = MAX( Element % BDOFs, Mesh % MaxBDOFs )
     Mesh % MaxNDOFs = MAX(Element % NDOFs / Element % TYPE % NumberOfNodes, &
         Mesh % MaxNDOFs)
   END DO

   DO i=1,Mesh % NumberOFBulkElements
     Element => Mesh % Elements(i)

     ! Set max element dofs here (because element size may have changed
     ! when edges and faces have been set). This is the absolute worst case.
     ! Element which has MaxElementDOFs may not even be present as a 
     ! real element
     Mesh % MaxElementDOFs = MAX( Mesh % MaxElementDOFs, &
          Element % TYPE % NumberOfNodes * Mesh % MaxNDOFs + &
          Element % TYPE % NumberOfEdges * Mesh % MaxEdgeDOFs + &
          Element % TYPE % NumberOfFaces * Mesh % MaxFaceDOFs + &
          Element % BDOFs, &
          Element % DGDOFs )

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


!------------------------------------------------------------------------------
  SUBROUTINE MeshStabParams( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    INTEGER :: i,n, istat
    LOGICAL :: stat, Stabilize, UseLongEdge
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

    CALL Info('MeshStabParams','Computing stabilization parameters',Level=7)
    CALL ResetTimer('MeshStabParams')

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('MeshStabParams','Mesh not associated')
    END IF
    
    IF ( Mesh % NumberOfNodes <= 0 ) RETURN

    Stabilize = .FALSE.
    
    DO i=1,CurrentModel % NumberOfSolvers
      Solver => CurrentModel % Solvers(i)
      IF ( ASSOCIATED( Mesh, Solver % Mesh ) ) THEN
        Stabilize = Stabilize .OR. &
            ListGetLogical( Solver % Values, 'Stabilize', Stat )
        Stabilize = Stabilize .OR. ListGetString( Solver % Values,  &
            'Stabilization Method', Stat )=='vms'
        Stabilize = Stabilize .OR.  ListGetString( Solver % Values, &
            'Stabilization Method', Stat )=='stabilized'
      END IF
    END DO

    Mesh % Stabilize = Stabilize 
    
    IF( ListGetLogical(CurrentModel % Simulation, &
        "Skip Mesh Stabilization",Stat) ) RETURN
    
    !IF( .NOT. Stabilize ) THEN
    !  CALL Info('MeshStabParams','No need to compute stabilization parameters',Level=10)      
    !  RETURN      
    !END IF
    
    CALL AllocateVector( Nodes % x, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % y, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % z, Mesh % MaxElementNodes )

    UseLongEdge = ListGetLogical(CurrentModel % Simulation, &
         "Stabilization Use Longest Element Edge",Stat)

    DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       n = Element % TYPE % NumberOfNodes
       Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
       Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
       Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
       IF ( Mesh % Stabilize ) THEN
          CALL StabParam( Element, Nodes,n, &
              Element % StabilizationMK, Element % hK, UseLongEdge=UseLongEdge)
       ELSE
          Element % hK = ElementDiameter( Element, Nodes, UseLongEdge=UseLongEdge)
       END IF
    END DO
 
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

    CALL CheckTimer('MeshStabParams',Level=7,Delete=.TRUE.)
!----------------------------------------------------------------------------
  END SUBROUTINE MeshStabParams
!------------------------------------------------------------------------------






!------------------------------------------------------------------------------
!> The quadratic mesh should be such that the center nodes lie roughly between
!> the corner nodes. This routine checks that this is actually the case.
!> The intended use for the routine is different kind of mesh related debugging.
!------------------------------------------------------------------------------
  SUBROUTINE InspectQuadraticMesh( Mesh, EnforceToCenter ) 
    
    TYPE(Mesh_t), TARGET :: Mesh
    LOGICAL, OPTIONAL :: EnforceToCenter

    LOGICAL :: Enforce
    INTEGER :: i,n,k,k1,k2,k3,ElemCode,ElemFamily,ElemDegree,ErrCount,TotCount
    REAL(KIND=dp) :: Center(3),Ref(3),Dist,Length
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: CenterMap(:,:)
    INTEGER, TARGET  :: TriangleCenterMap(3,3), QuadCenterMap(4,3), &
        TetraCenterMap(6,3), BrickCenterMap(12,3), WedgeCenterMap(9,3), PyramidCenterMap(8,3) 
    
    CALL Info('InspectQuadraticMesh','Inspecting quadratic mesh for outliers')
    CALL Info('InspectQuadraticMesh','Number of nodes: '//I2S(Mesh % NumberOfNodes),Level=8)
    CALL Info('InspectQuadraticMesh','Number of bulk elements: '&
        //I2S(Mesh % NumberOfBulkElements),Level=8)
    CALL Info('InspectQuadraticMesh','Number of boundary elements: '&
        //I2S(Mesh % NumberOfBoundaryElements),Level=8)


    IF( PRESENT( EnforceToCenter ) ) THEN
      Enforce = EnforceToCenter
    ELSE
      Enforce = .FALSE.
    END IF

    TriangleCenterMap(1,:) = [ 1, 2, 4]
    TriangleCenterMap(2,:) = [ 2, 3, 5]
    TriangleCenterMap(3,:) = [ 3, 1, 6]
    
    QuadCenterMap(1,:) = [ 1, 2, 5]
    QuadCenterMap(2,:) = [ 2, 3, 6]
    QuadCenterMap(3,:) = [ 3, 4, 7]
    QuadCenterMap(4,:) = [ 4, 1, 8]
    
    TetraCenterMap(1,:) = [ 1, 2, 5]
    TetraCenterMap(2,:) = [ 2, 3, 6]
    TetraCenterMap(3,:) = [ 3, 1, 7]
    TetraCenterMap(4,:) = [ 1, 4, 8]
    TetraCenterMap(5,:) = [ 2, 4, 9]
    TetraCenterMap(6,:) = [ 3, 4, 10]

    BrickCenterMap(1,:) = [ 1, 2,  9 ]
    BrickCenterMap(2,:) = [ 2, 3,  10 ]
    BrickCenterMap(3,:) = [ 3, 4,  11 ]
    BrickCenterMap(4,:) = [ 4, 1,  12 ]
    BrickCenterMap(5,:) = [ 1, 5,  13 ]
    BrickCenterMap(6,:) = [ 2, 6,  14 ]
    BrickCenterMap(7,:) = [ 3, 7,  15 ]
    BrickCenterMap(8,:) = [ 4, 8,  16 ]
    BrickCenterMap(9,:) = [ 5, 6,  17 ]
    BrickCenterMap(10,:) = [ 6, 7, 18 ]
    BrickCenterMap(11,:) = [ 7, 8, 19 ]
    BrickCenterMap(12,:) = [ 8, 5, 20 ]
    
    WedgeCenterMap(1,:) = [ 1, 2, 7 ]
    WedgeCenterMap(2,:) = [ 2, 3, 8 ]
    WedgeCenterMap(3,:) = [ 3, 1, 9 ]
    WedgeCenterMap(4,:) = [ 4, 5, 10 ]
    WedgeCenterMap(5,:) = [ 5, 6, 11 ]
    WedgeCenterMap(6,:) = [ 6, 4, 12 ]
    WedgeCenterMap(7,:) = [ 1, 4, 13 ]
    WedgeCenterMap(8,:) = [ 2, 5, 14 ]
    WedgeCenterMap(9,:) = [ 3, 6, 15 ]
    
    PyramidCenterMap(1,:) = [ 1,2,6 ]
    PyramidCenterMap(2,:) = [ 2,3,7 ]
    PyramidCenterMap(3,:) = [ 3,4,8 ]
    PyramidCenterMap(4,:) = [ 4,1,9 ]
    PyramidCenterMap(5,:) = [ 1,5,10 ]
    PyramidCenterMap(6,:) = [ 2,5,11 ]
    PyramidCenterMap(7,:) = [ 3,5,12 ]
    PyramidCenterMap(8,:) = [ 4,5,13 ]
    
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z
    
    !   Loop over elements:
    !   -------------------
    ErrCount = 0
    TotCount = 0

    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(i)

      ElemCode = Element % TYPE % ElementCode 
      ElemFamily = ElemCode / 100
      ElemDegree = Element % TYPE % BasisFunctionDegree
      
      ! Only check quadratic elements!
      IF( ElemDegree /= 2 ) CYCLE
      
      SELECT CASE( ElemFamily ) 

      CASE(3)
        n = 3
        CenterMap => TriangleCenterMap
        
      CASE(4)
        n = 4
        CenterMap => QuadCenterMap
        
      CASE(5)
        n = 6
        CenterMap => TetraCenterMap
        
      CASE(6)
        n = 8
        CenterMap => PyramidCenterMap
        
      CASE(7)
        n = 9
        CenterMap => WedgeCenterMap
        
      CASE(8)
        n = 12
        CenterMap => BrickCenterMap
        
      CASE DEFAULT
        CALL Fatal('InspectQuadraticMesh','Element type '//I2S(ElemCode)//' not implemented!')

      END SELECT
      
      !      Loop over every edge of every element:
      !      --------------------------------------
       DO k=1,n
         k1 = Element % NodeIndexes( CenterMap(k,1) )
         k2 = Element % NodeIndexes( CenterMap(k,2) )
         k3 = Element % NodeIndexes( CenterMap(k,3) )
         
         Center(1) = ( x(k1) + x(k2) ) / 2.0_dp
         Center(2) = ( y(k1) + y(k2) ) / 2.0_dp
         Center(3) = ( z(k1) + z(k2) ) / 2.0_dp

         Ref(1) = x(k3)
         Ref(2) = y(k3) 
         Ref(3) = z(k3)

         Length = SQRT( (x(k1) - x(k2))**2.0 + (y(k1) - y(k2))**2.0 + (z(k1) - z(k2))**2.0 )
         Dist = SQRT( SUM( (Center - Ref)**2.0 ) )

         TotCount = TotCount + 1
         IF( Dist > 0.01 * Length ) THEN
           ErrCount = ErrCount + 1
           PRINT *,'Center Displacement:',i,ElemCode,n,k,Dist/Length
         END IF

         IF( Enforce ) THEN
           x(k3) = Center(1)
           y(k3) = Center(2)
           z(k3) = Center(3)
         END IF

       END DO
     END DO
         
     IF( TotCount > 0 ) THEN
       CALL Info('InspectQuadraticMesh','Number of outlier nodes is '&
           //I2S(ErrCount)//' out of '//I2S(TotCount),Level=6)
     ELSE
       CALL Info('InspectQuadraticMesh','No quadratic elements to inspect',Level=8)
     END IF

  END SUBROUTINE InspectQuadraticMesh





  !---------------------------------------------------------------------------
  SUBROUTINE PolynomBoundaryFit(Mesh, PParams, BCind, Ndeg, FitParams, PatchHeight )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    TYPE(Valuelist_t), POINTER :: PParams
    INTEGER, OPTIONAL :: BCind
    INTEGER :: Ndeg
    REAL(KIND=dp) :: FitParams(:)
    REAL(KIND=dp), POINTER :: PatchHeight(:)

    REAL(KIND=dp), POINTER :: rArray(:,:)
    LOGICAL :: Found
    REAL(KIND=dp), POINTER :: pArray(:,:)    
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    LOGICAL, ALLOCATABLE :: ActiveNode(:)
    REAL(KIND=dp), ALLOCATABLE :: AngleSum(:), pHeight(:), Weights(:)
    TYPE(Nodes_t) :: Nodes
    INTEGER :: n,nd,dim,Corners(4)
    
    pArray => ListGetConstRealArray( PParams,'Patch Height Basis',Found )

    IF(.NOT. Found ) THEN
      dim = 3
      x => Mesh % Nodes % x
      y => Mesh % Nodes % y
      z => Mesh % Nodes % z

      n = Mesh % NumberOfNodes
      ALLOCATE(ActiveNode(n), AngleSum(n), Weights(n) )
      ActiveNode = .FALSE.
      AngleSum = 0.0_dp
      Weights = 0.0_dp
      
      CALL FindBoundaryCorners()      
      
      CALL SetBoundaryWeights()
      
      CALL FitBoundaryPatch()

      NULLIFY(rArray)
      ALLOCATE(rArray(nd,1))
      rArray = 0.0_dp

      rArray(1:4,1) = Nodes % x(1:4)
      CALL ListAddConstRealArray( PParams,'Patch Corners x',4,1,rArray) 
      rArray(1:4,1) = Nodes % y(1:4)
      CALL ListAddConstRealArray( PParams,'Patch Corners y',4,1,rArray)
      rArray(1:4,1) = Nodes % z(1:4)
      CALL ListAddConstRealArray( PParams,'Patch Corners z',4,1,rArray)
      rArray(1:nd,1) = pheight(1:nd)
      CALL ListAddConstRealArray( PParams,'Patch Height Basis',nd,1,rArray)
      DEALLOCATE(pheight)
      
      pArray => ListGetConstRealArray( PParams,'Patch Height Basis',Found )
    END IF

    ALLOCATE(PatchHeight(SIZE(pArray,1)))
    PatchHeight = pArray(:,1)    
    pArray => ListGetConstRealArray( PParams,'Patch Corners x',UnfoundFatal=.TRUE. )    
    FitParams(1:4) = pArray(1:4,1)
    pArray => ListGetConstRealArray( PParams,'Patch Corners y',UnfoundFatal=.TRUE. )
    FitParams(5:8) = pArray(1:4,1)
    pArray => ListGetConstRealArray( PParams,'Patch Corners z',UnfoundFatal=.TRUE. )
    FitParams(9:12) = pArray(1:4,1)
      
  CONTAINS
    
    ! Found the four courners of the patch. It is assumed that they are the ones
    ! with the smallest angle. Typically that would be 90 degs.
    !---------------------------------------------------------------------------------
    SUBROUTINE FindBoundaryCorners()

      INTEGER :: t,t1,t2,i,j,k,i1,i2,j1,j2
      REAL(KIND=dp) :: v1(3),v2(3),phi,Angles(4),dist,maxdist
      TYPE(Element_t), POINTER :: Element
      
      t1 = Mesh % NumberOfBulkElements
      t2 = Mesh % NumberOfBoundaryElements 
      
      DO t=t1+1,t1+t2
        Element => Mesh % Elements(t)
        IF ( Element % BoundaryInfo % Constraint /= CurrentModel % BCs(BCind) % Tag ) CYCLE
        
        n  = MODULO(Element % TYPE % ElementCode, 100)
        IF(n < 3 .OR. n > 4 ) THEN
          CALL Fatal('PolynomBoundaryFit','2D polynom can only bet fitted on 2D elements!')
        END IF
        
        DO i=1,n
          i1 = MODULO(i-2,n)+1
          i2 = MODULO(i,n)+1
          j = Element % NodeIndexes(i)
          j1 = Element % NodeIndexes(i1)
          j2 = Element % NodeIndexes(i2)
          v1(1) = x(j1)-x(j)
          v1(2) = y(j1)-y(j)
          v1(3) = z(j1)-z(j)
          v2(1) = x(j2)-x(j)
          v2(2) = y(j2)-y(j)
          v2(3) = z(j2)-z(j)
          v1 = v1 / SQRT(SUM(v1**2))
          v2 = v2 / SQRT(SUM(v2**2))        
          phi = ACOS(SUM(v1*v2))        
          AngleSum(j) = AngleSum(j) + phi
          ActiveNode(j) = .TRUE.
        END DO
      END DO

      Angles = HUGE(phi)
      DO j=1,4
        k = MINLOC(AngleSum, dim = 1, Mask = ActiveNode )
        Corners(j) = k
        Angles(j) = AngleSum(k) 
        ! Eliminate the minimum angle and repeat to find the next smallest angle. 
        AngleSum(k) = 3*PI
      END DO

      IF( InfoActive(20 ) ) THEN
        Angles = ( 180_dp / PI ) * Angles
        PRINT *,'Patch element corners:',Corners
        PRINT *,'Patch element angles:',Angles
      END IF

      ! Find the two nodes furthers apart.
      maxdist = 0.0_dp
      DO i=1,4
        DO j=i+1,4
          v1(1) = x(Corners(j))-x(Corners(i))
          v1(2) = y(Corners(j))-y(Corners(i))
          v1(3) = z(Corners(j))-z(Corners(i))
          dist = SQRT(SUM(v1**2))
          IF(dist > maxdist) THEN
            i1 = i
            i2 = j
            maxdist = dist
          END IF
        END DO
      END DO

      ! Swap the nodes furthest apart so that they are nodes 1 and 3 (always: i2>i1)
      ! Then (1-2) and (1-4) create two basis vectors for the plane.
      IF(i1 /= 1 ) THEN
        k = Corners(1)
        Corners(1) = Corners(i1)
        Corners(i1) = k
      END IF        
      IF(i2 /= 3) THEN
        k = Corners(3)
        Corners(3) = Corners(i2)
        Corners(i2) = k
      END IF
      
    END SUBROUTINE FindBoundaryCorners   


    ! We want to set the value at nodes, not at integration points. However, we need to sum 
    ! up the integration weights to the nodes.
    !---------------------------------------------------------------------------------------
    SUBROUTINE SetBoundaryWeights()
      TYPE(Element_t), POINTER :: sElement
      TYPE(Nodes_t) :: sNodes
      INTEGER :: t,t1,t2,i,n
      INTEGER, POINTER :: sIndexes(:)
      REAL(KIND=dp) :: Basis(4), detJ
      TYPE(GaussIntegrationPoints_t) :: IP
      LOGICAL :: stat
            
      t1 = Mesh % NumberOfBulkElements
      t2 = Mesh % NumberOfBoundaryElements 
      
      DO t=t1+1,t1+t2
        sElement => Mesh % Elements(t)        
        IF ( sElement % BoundaryInfo % Constraint /= CurrentModel % BCs(BCind) % Tag ) CYCLE

        sIndexes => sElement % NodeIndexes
        n  = sElement % TYPE % NumberOfNodes

        IP = GaussPoints( sElement )
        CALL CopyElementNodesFromMesh( sNodes, Mesh, n, sIndexes)
        
        DO i=1,IP % n
          stat = ElementInfo( sElement, sNodes, IP % U(i), IP % V(i), &
              IP % W(i), detJ, Basis )
          Weights(sIndexes) = Weights(sIndexes) + IP % s(i) * detJ * Basis(1:n)
        END DO
      END DO

      IF( InfoActive(20) ) THEN
        PRINT *,'Sum of Weights on element patch:',SUM(Weights)
      END IF
        
    END SUBROUTINE SetBoundaryWeights


    ! Fit patch to the height data given on each node of the patch.
    !--------------------------------------------------------------
    SUBROUTINE FitBoundaryPatch()
      TYPE(Element_t), TARGET :: Element
      TYPE(Element_t), POINTER :: pElement
      INTEGER :: n,i,j,k,q,np,edofs
      REAL(KIND=dp) :: c1(3), c2(3), c4(3), normal(3), v1(3), v2(3), &
          u,v,w, weight, detJ, norm_proj, dir(4)
      REAL(KIND=dp), ALLOCATABLE :: MASS(:,:), FORCE(:), Basis(:)
      LOGICAL :: Erroneous, Stat, Invert, Serendipity
      INTEGER :: pivot(50)
      TYPE(GaussIntegrationPoints_t) :: IP

      ! Define parameters for p-element patch. 
      n = 4
      np = (ndeg+1)**2      
      edofs = ndeg - 1      
      nd = n*(1+edofs)
      Serendipity = .TRUE.
      
      IF(.NOT. ASSOCIATED(Nodes % x) ) THEN
        ALLOCATE(Nodes % x(nd), Nodes % y(nd), Nodes % z(nd), Basis(nd)) 
        Nodes % x = 0.0_dp; Nodes % y = 0.0_dp; Nodes % z = 0.0_dp
        Basis = 0.0_dp
      END IF        

      DO i=1,n        
        Nodes % x(i) = x(Corners(i))
        Nodes % y(i) = y(Corners(i))
        Nodes % z(i) = z(Corners(i))
      END DO
      
      ! Creat basis vectors for the element assuming that it can be in a plane.  
      c1(1) = x(Corners(1))
      c1(2) = y(Corners(1))
      c1(3) = z(Corners(1))
      c2(1) = x(Corners(2))
      c2(2) = y(Corners(2))
      c2(3) = z(Corners(2))
      c4(1) = x(Corners(4))
      c4(2) = y(Corners(4))
      c4(3) = z(Corners(4))
      
      normal = NormalDirection(c2-c1,c4-c1)
            
      Element % TYPE => GetElementType(404)
      pElement => Element
      
      ALLOCATE(MASS(nd,nd),FORCE(nd),pheight(nd))
      MASS = 0.0_dp
      FORCE = 0.0_dp
      pheight = 0.0_dp
      Weight = 1.0_dp
      
      IP = GaussPoints( pElement, np = np, PReferenceElement = .TRUE.)

      ! Currently equal weight for all nodes.
      DO i=1,Mesh % NumberOfNodes
        IF(.NOT. ActiveNode(i)) CYCLE

        v1(1) = x(i) 
        v1(2) = y(i) 
        v1(3) = z(i) 
        
        norm_proj = SUM((v1-c1)*normal)
        
        ! We can only find the integration points on the plane defined by the superelement. 
        v2 = v1 - norm_proj * normal

        CALL GlobalToLocal( u,v,w,v2(1),v2(2),v2(3),pElement,Nodes )

        weight = weights(i)

        ! This is minimal quadrilateral p-element on-the-fly without any excess definions needed
        q = n
        CALL QuadNodalPBasisAll(u, v, basis) 
        DO j=1,4
          invert = (j==4)
          DO k=1,edofs
            q = q + 1
            ! Get values of basis functions for edge=j and j=k+1 by parity
            IF (Serendipity) THEN
              Basis(q) = SD_QuadEdgePBasis(j,k+1,u,v,invert)
            ELSE
              Basis(q) = QuadEdgePBasis(j,k+1,u,v,invert)
            END IF
          END DO
        END DO

        ! Create equation involving mass matrix that solves for the coordinates at the p-dofs
        DO q=1,nd
          MASS(1:nd,q) = MASS(1:nd,q) + Weight * Basis(1:nd) * Basis(q) 
        END DO        
        FORCE(1:nd) = FORCE(1:nd) + Weight * Basis(1:nd) * norm_proj
        
        DO j=1,4
          IF(Corners(j) == i) dir(j) = norm_proj
        END DO        
      END DO

      ! Set dirichlet conditions for the corners
      DO j=1,4
        MASS(j,1:nd) = 0.0_dp
        MASS(j,j) = 1.0_dp
        FORCE(j) = dir(j)
      END DO
      
      CALL LUdecomp(MASS,nd,pivot,Erroneous)
      IF (Erroneous) CALL Fatal('FitBoundaryPatch', 'LU-decomposition fails')      
      pheight = FORCE
      CALL LUSolve(nd,MASS,pheight,pivot)

      DEALLOCATE(MASS,FORCE)
      
    END SUBROUTINE FitBoundaryPatch
          
  END SUBROUTINE PolynomBoundaryFit

  
  
  SUBROUTINE FollowCurvedBoundary(Model, Mesh, SetP )
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), TARGET :: Mesh 
    LOGICAL :: SetP

    LOGICAL :: Found
    REAL(KIND=dp) :: FitParams(12)
    REAL(KIND=dp), POINTER :: normheight(:)
    INTEGER :: Mode, bc_ind, dim, ndeg
    TYPE(ValueList_t), POINTER :: BC

    IF(.NOT. ListCheckPrefixAnyBC( Model,'Follow') ) RETURN

    dim = Mesh % MeshDim
    
    FitParams = 0
    DO bc_ind = 1, Model % NumberOfBCs
      BC => Model % BCs(bc_ind) % Values
      IF( ListGetLogical(BC,'Follow Circle Boundary', Found ) ) THEN
        CALL CylinderFit(Mesh, BC, bc_ind, 2, FitParams ) 
        Mode = 1        
      ELSE IF( ListGetLogical(BC,'Follow Cylinder Boundary', Found ) ) THEN
        CALL CylinderFit(Mesh, BC, bc_ind, dim, FitParams) 
        Mode = 2        
      ELSE IF( ListGetLogical(BC,'Follow Sphere Boundary', Found ) ) THEN
        CALL SphereFit(Mesh, BC, bc_ind, FitParams ) 
        Mode = 3        
      ELSE IF( ListGetLogical(BC,'Follow Function Boundary', Found ) ) THEN
        IF(.NOT. ListCheckPresent(BC,'Surface Function') ) THEN
          CALL Fatal('FollowCurvedBoundary','We need "Surface Function" to follow!')
        END IF
        Mode = 4        
      ELSE IF( ListGetLogical(BC,'Follow Toroid Boundary', Found ) ) THEN
        CALL TorusFit(Mesh, BC, bc_ind, FitParams ) 
        Mode = 5        
      ELSE IF( ListCheckPresent(BC,'Follow Polynom Boundary' ) ) THEN
        ndeg = ListGetInteger( BC,'Follow Polynom Boundary', Found ) 
        CALL PolynomBoundaryFit(Mesh, BC, bc_ind, Ndeg, FitParams, normheight ) 
        Mode = 6        
      ELSE
        Mode = 0
      END IF
      
      IF(Mode > 0 ) THEN
        CALL Info('FollowCurvedBoundary','Setting BC '//I2S(bc_ind)//&
            ' to follow curved boundary in mode '//I2S(Mode),Level=7)
        CALL SetCurvedBoundary()
      END IF
    END DO

    
  CONTAINS

    
    ! We have fitted a p-element patch to a rectangular boundary.
    ! Now apply if to each element of the boundary.
    !--------------------------------------------------------------
    FUNCTION PatchElementApply(v1) RESULT( v2 )
      REAL(KIND=dp) :: v1(3)
      REAL(KIND=dp) :: v2(3) 

      REAL(KIND=dp) :: c1(3),c2(3),c4(3),normal(3),norm_proj,u,v,w,h
      INTEGER :: n,np,nd,q,edofs,i,j,k
      TYPE(Element_t), TARGET :: Element
      TYPE(Element_t), POINTER :: pElement
      TYPE(Nodes_t), SAVE :: Nodes
      REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:)
      LOGICAL :: Serendipity, invert
            
      ! Create basis functions using the corners
      c1 = FitParams([1,5,9])
      c2 = FitParams([2,6,10])
      c4 = FitParams([4,8,12])
      
      ! Remove normal components so we are in plane
      ! We can only find the integration points on the plane defined by the superelement. 
      normal = NormalDirection(c2-c1,c4-c1)
      norm_proj = SUM((v1-c1)*normal)
      v2 = v1 - norm_proj * normal

      ! Parameters of the p-element
      n = 4
      np = (ndeg+1)**2      
      edofs = ndeg - 1      
      nd = n*(1+edofs)
      Serendipity = .TRUE.

      IF(.NOT. ASSOCIATED(Nodes % x)) THEN
        ALLOCATE(Nodes % x(n), Nodes % y(n), Nodes % z(n),Basis(nd))
      END IF
      Nodes % x = FitParams(1:4)
      Nodes % y = FitParams(5:8)
      Nodes % z = FitParams(9:12)
      Basis = 0.0_dp

      ! Find local coordinates of the node it the patch element. 
      Element % TYPE => GetElementType(404)
      pElement => Element

      ! Give the global coordinates in loca coordinates of the patch element. 
      CALL GlobalToLocal( u,v,w,v2(1),v2(2),v2(3),pElement,Nodes )
      
      ! This is minimal quadrilateral p-element on-the-fly without any excess definions needed
      ! Given the local coordinates find the basis function values at the point.
      q = n
      CALL QuadNodalPBasisAll(u, v, basis) 
      DO j=1,4
        invert = (j==4)
        DO k=1,edofs
          q = q + 1
          ! Get values of basis functions for edge=j and j=k+1 by parity
          IF (Serendipity) THEN
            Basis(q) = SD_QuadEdgePBasis(j,k+1,u,v,invert)
          ELSE
            Basis(q) = QuadEdgePBasis(j,k+1,u,v,invert)
          END IF
        END DO
      END DO

      ! Get the updated height and return the new coordinates. 
      h = SUM( Basis(1:nd) * normheight(1:nd) )
      v2 = v2 + h * Normal
      
    END FUNCTION PatchElementApply

    
          
!------------------------------------------------------------------------------
    SUBROUTINE SetCurvedBoundary()
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: R, Rminor, r1, rat, f, gradf(3)
      REAL(KIND=dp) :: Nrm(3), Tngt1(3), Tngt2(3), Orig(3), Coord(3), NtCoord(3), PlaneCoord(3)
      INTEGER :: i,j,k,l,t,n,elem
      LOGICAL, POINTER :: DoneNode(:)
      TYPE(Element_t), POINTER :: Element
      LOGICAL :: Parallel 
      TYPE(ParallelInfo_t), POINTER :: ParallelInfo
      
      IF( Mode == 1 ) THEN  ! circle
        Orig(1:2) = FitParams(1:2)
        Orig(3) = 0.0_dp
        R = FitParams(3)
        IF( InfoActive(25) .AND. ParEnv % MyPe == 0) PRINT *,'Circle Params:',FitParams(1:3)                        
      ELSE IF( Mode == 2 ) THEN  ! cylinder 
        Orig(1:3) = FitParams(1:3)
        Nrm(1:3) = FitParams(4:6)        
        R = FitParams(7)
        IF( InfoActive(25) .AND. ParEnv % MyPe == 0) PRINT *,'Cylinder Params:',FitParams(1:7)        
        CALL TangentDirections(Nrm, Tngt1, Tngt2 ) 
      ELSE IF( Mode == 3 ) THEN ! sphere
        Orig(1:3) = FitParams(1:3)
        Nrm = 0.0_dp
        R = FitParams(4)
        IF( InfoActive(25) .AND. ParEnv % MyPe == 0) PRINT *,'Sphere Params:',FitParams(1:4)                                
      ELSE IF( Mode == 4 ) THEN
        Orig = 0.0_dp        
      ELSE IF( Mode == 5 ) THEN  ! torus
        Orig(1:3) = FitParams(1:3)
        Nrm(1:3) = FitParams(4:6)        
        R = FitParams(7)
        Rminor = FitParams(8)
        IF( InfoActive(25) .AND. ParEnv % MyPe == 0) PRINT *,'Torus Params:',FitParams(1:8)        
        CALL TangentDirections(Nrm, Tngt1, Tngt2 ) 
      ELSE IF( Mode == 6 ) THEN
        Orig = 0.0_dp
      END IF
      
      Parallel = ( ParEnv % PEs > 1 .AND. .NOT. Mesh % SingleMesh )

      PRINT *,'SetP:',SetP
      
      IF(.NOT. SetP) THEN
        ALLOCATE( DoneNode(Mesh % NumberOfNodes))
        DoneNode = .FALSE.
        
        DO elem=Mesh % NumberOfBulkElements+1, &
            Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(elem)
          IF ( Element % BoundaryInfo % Constraint &
              == Model % BCs(bc_ind) % Tag ) THEN      
            n = Element % TYPE % NumberOfNodes          
            DoneNode(Element % NodeIndexes(1:n)) = .TRUE.
          END IF
        END DO

        IF( Parallel ) THEN
          ParallelInfo => Mesh % ParallelInfo 
          CALL CommunicateParallelSystemTag(ParallelInfo,Ltag = DoneNode)
        END IF

        DO j=1, Mesh % NumberOfNodes
          IF( .NOT. DoneNode(j) ) CYCLE

          Coord(1) = Mesh % Nodes % x(j) - Orig(1)           
          Coord(2) = Mesh % Nodes % y(j) - Orig(2)
          Coord(3) = Mesh % Nodes % z(j) - Orig(3)
          
          SELECT CASE( Mode )
          CASE( 1 ) ! circle 
            rat = R / SQRT(SUM(Coord(1:2)**2))
            Coord(1:2) = rat*Coord(1:2)
          CASE( 2 ) ! cylinder
            NtCoord(1) = SUM(Nrm*Coord)
            NtCoord(2) = SUM(Tngt1*Coord)
            NtCoord(3) = SUM(Tngt2*Coord)
            rat = R / SQRT(SUM(NtCoord(2:3)**2))
            NtCoord(2:3) = rat*NtCoord(2:3)
            Coord = NtCoord(1)*Nrm + NtCoord(2)*Tngt1 + NtCoord(3)*Tngt2
          CASE( 3 ) ! sphere 
            rat = R / SQRT(SUM(Coord(1:3)**2))
            Coord(1:3) = rat*Coord(1:3)
          CASE( 4 ) ! analytical function
            ! For now we fix Newton's iteration to three...
            DO i=1,3
              f = ListGetFunVec( BC,'Surface Function', Coord(1:dim), dim, DfDx=gradf(1:dim) )
              Coord(1:dim) = Coord(1:dim) - f*gradf(1:dim)/(SUM(gradf(1:dim)**2))
            END DO
          CASE( 5 ) ! torus
            NtCoord(1) = SUM(Nrm*Coord)
            NtCoord(2) = SUM(Tngt1*Coord)
            NtCoord(3) = SUM(Tngt2*Coord)

            PlaneCoord(1) = 0.0_dp
            r1 = SQRT(SUM(NtCoord(2:3)**2))
            PlaneCoord(2:3) = NtCoord(2:3)*R/r1

            rat = Rminor / SQRT((r1-R)**2 + NtCoord(1)**2)            
            NtCoord = rat * (NtCoord-PlaneCoord) + PlaneCoord

            Coord = NtCoord(1)*Nrm + NtCoord(2)*Tngt1 + NtCoord(3)*Tngt2

          CASE( 6 ) 
            Coord = PatchElementApply(Coord)
            
          END SELECT
          
          Mesh % Nodes % x(j) = Coord(1) + Orig(1)
          Mesh % Nodes % y(j) = Coord(2) + Orig(2)
          Mesh % Nodes % z(j) = Coord(3) + Orig(3)
        END DO
        DEALLOCATE(DoneNode)
      END IF
        
      IF( SetP ) THEN
        DO elem=Mesh % NumberOfBulkElements+1, &
            Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(elem)
          IF ( Element % BoundaryInfo % Constraint &
              /= Model % BCs(bc_ind) % Tag ) CYCLE          
          n = Element % TYPE % NumberOfNodes
          
          BLOCK 
            REAL(KIND=dp) :: Weight
            REAL(KIND=dp) :: Basis(50),DetJ
            REAL(KIND=dp) :: MASS(50,50), FORCE(3,50), x(50), Coord0(3)
            LOGICAL :: Stat, Erroneous
            INTEGER :: nd,i,t,p,q
            INTEGER, TARGET :: Indexes(50)
            INTEGER :: pivot(50)
            INTEGER, POINTER :: pIndexes(:)
            TYPE(GaussIntegrationPoints_t) :: IP
            TYPE(Nodes_t), SAVE :: Nodes

            pIndexes => Indexes 
            Nd = mGetElementDOFs( pIndexes, Element, CurrentModel % Solver )          
                        
            ! Only if we have really p-elements is there a need to consider the curved shape
            IF(Nd == n ) CYCLE

            CALL CopyElementNodesFromMesh( Nodes, Mesh, n, pIndexes)

            MASS = 0._dp
            FORCE = 0._dp

            IP = GaussPoints( Element )
            
            DO t=1,IP % n
              stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                  IP % W(t), detJ, Basis )
              Weight = IP % s(t) * DetJ

              ! Current nodal value at integration point does not consider p-dofs
              Coord(1) = SUM( Nodes % x(1:n) * Basis(1:n) )
              Coord(2) = SUM( Nodes % y(1:n) * Basis(1:n) )
              Coord(3) = SUM( Nodes % z(1:n) * Basis(1:n) )
              Coord0 = Coord

              Coord = Coord - Orig
              SELECT CASE( Mode )
              CASE( 1 ) 
                rat = R / SQRT(SUM(Coord(1:2)**2))
                Coord(1:2) = rat * Coord(1:2)
              CASE( 2 )
                ! Local coordinates in nt-system
                NtCoord(1) = SUM(Nrm*Coord)
                NtCoord(2) = SUM(Tngt1*Coord)
                NtCoord(3) = SUM(Tngt2*Coord)
                ! Ratio between current and desired radius
                rat = R / SQRT(SUM(NtCoord(2:3)**2))
                NtCoord(2:3) = rat * NtCoord(2:3)
                Coord = NtCoord(1)*Nrm + NtCoord(2)*Tngt1 + NtCoord(3)*Tngt2
              CASE( 3 ) 
                rat = R / SQRT(SUM(Coord(1:3)**2))
                Coord(1:3) = rat * Coord(1:3)
              CASE( 4 ) 
                DO i=1,3
                  f = ListGetFunVec( BC,'Surface Function', Coord(1:dim), dim, DfDx=gradf(1:dim) )
                  Coord(1:dim) = Coord(1:dim) - f*gradf(1:dim)/(SUM(gradf(1:dim)**2))            
                END DO
              CASE( 5 ) ! torus
                NtCoord(1) = SUM(Nrm*Coord)
                NtCoord(2) = SUM(Tngt1*Coord)
                NtCoord(3) = SUM(Tngt2*Coord)

                PlaneCoord(1) = 0.0_dp
                r1 = SQRT(SUM(NtCoord(2:3)**2))
                PlaneCoord(2:3) = NtCoord(2:3)*R/r1

                rat = Rminor / SQRT((r1-R)**2 + NtCoord(1)**2)            
                NtCoord = rat * (NtCoord-PlaneCoord) + PlaneCoord
                Coord = NtCoord(1)*Nrm + NtCoord(2)*Tngt1 + NtCoord(3)*Tngt2
              CASE( 6 ) 
                Coord = PatchElementApply(Coord)
              END SELECT

              Coord = Coord + Orig
              ! Solve for desired coordinate displacement rather than absolute coordinate value
              Coord = Coord - Coord0
              
              ! Create equation involving mass matrix that solves for the coordinates at the p-dofs
              DO q=1,nd
                MASS(1:nd,q) = MASS(1:nd,q) + Weight * Basis(1:nd) * Basis(q) 
              END DO

              DO i=1,dim
                FORCE(i,1:nd) = FORCE(i,1:nd) + Weight * Basis(1:nd) * Coord(i) 
              END DO
            END DO

            ! Set Dirichlet conditions for the nodal coordinate displacements
            DO i=1,n
              MASS(i,1:nd) = 0.0_dp
              MASS(i,i) = 1.0_dp
              FORCE(:,i) = 0.0_dp
            END DO
            
            CALL LUdecomp(MASS,nd,pivot,Erroneous)
            IF (Erroneous) THEN
              PRINT *,'Element:',elem,ip % n,nd,n,dim
              DO i=1,dim
                PRINT *,'FORCE',i,FORCE(i,1:nd)
              END DO
              DO i=n+1,nd
                PRINT *,'MASS',i,MASS(i,1:nd)
              END DO
              CALL Fatal('SetCurvedBoundary', 'LU-decomposition fails')
            END IF
              
            DO i=1,dim          
              x(1:nd) = FORCE(i,1:nd)
              CALL LUSolve(nd,MASS,x,pivot)
              
              SELECT CASE(i)
              CASE(1)
                Mesh % Nodes % x(Indexes(n+1:nd)) = x(n+1:nd) 
              CASE(2)
                Mesh % Nodes % y(Indexes(n+1:nd)) = x(n+1:nd) 
              CASE(3)
                Mesh % Nodes % z(Indexes(n+1:nd)) = x(n+1:nd) 
              END SELECT
            END DO
            
          END BLOCK
        END DO
      END IF
        
    END SUBROUTINE SetCurvedBoundary
!------------------------------------------------------------------------------
  END SUBROUTINE FollowCurvedBoundary

  
  
  !------------------------------------------------------------------------------------------------
  !> Finds nodes for which CandNodes are True such that their mutual distance is somehow
  !> maximized. We first find lower left corner, then the node that is furthest apart from it,
  !> and continue as long as there are nodes to find. Typically we would be content with two nodes
  !> on a line, three nodes on a plane, and four nodes on a volume.
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE FindExtremumNodes(Mesh,CandNodes,NoExt,Inds)
    TYPE(Mesh_t), TARGET :: Mesh
    LOGICAL, ALLOCATABLE :: CandNodes(:)
    INTEGER :: NoExt
    INTEGER, TARGET :: Inds(:)

    REAL(KIND=dp) :: Coord(3),dCoord(3),dist,MinDist,MaxDist
    REAL(KIND=dp), ALLOCATABLE :: SetCoord(:,:)
    INTEGER :: i,j,k
    
    ALLOCATE( SetCoord(NoExt,3) )
    SetCoord = 0.0_dp
    Inds = 0
    
    ! First find the lower left corner
    MinDist = HUGE(MinDist) 
    DO i=1, Mesh % NumberOfNodes
      IF(.NOT. CandNodes(i) ) CYCLE
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
      Dist = SUM( Coord )
      IF( Dist < MinDist ) THEN
        Inds(1) = i
        MinDist = Dist
        SetCoord(1,:) = Coord
      END IF
    END DO
    
    ! Find more points such that their minimum distance to the previous point(s)
    ! is maximized.
    DO j=2,NoExt
      ! The maximum minimum distance of any node from the previously defined nodes
      MaxDist = 0.0_dp
      DO i=1, Mesh % NumberOfNodes
        IF(.NOT. CandNodes(i) ) CYCLE
        Coord(1) = Mesh % Nodes % x(i)
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)
        
        ! Minimum distance from the previously defined nodes
        MinDist = HUGE(MinDist)
        DO k=1,j-1
          dCoord = SetCoord(k,:) - Coord
          Dist = SUM( dCoord**2 )          
          MinDist = MIN( Dist, MinDist )
        END DO
        
        ! If the minimum distance is greater than in any other node, choose this
        IF( MaxDist < MinDist ) THEN
          MaxDist = MinDist 
          Inds(j) = i
          SetCoord(j,:) = Coord
        END IF
      END DO
    END DO

    IF( InfoActive(30) ) THEN
      PRINT *,'Extremum Inds:',Inds
      DO i=1,NoExt
        PRINT *,'Node:',Inds(i),SetCoord(i,:)
      END DO
    END IF
      
  END SUBROUTINE FindExtremumNodes
    

    
  ! This creates a projector that integrates over the BCs on the boundary such that
  ! an integral constraint may be applied on it. For example, we could set the
  ! incoming flow without actually setting the profile.
  !--------------------------------------------------------------------------------------
  FUNCTION IntegralProjector(Model, Mesh, BCInd, IsBodyForce ) RESULT ( Projector )

    TYPE(Model_t) :: Model  
    TYPE(Mesh_t), TARGET :: Mesh
    INTEGER :: BCInd
    LOGICAL :: IsBodyForce
    TYPE(Matrix_t), POINTER :: Projector
        
    REAL(KIND=dp) :: area
    TYPE(ValueList_t), POINTER :: BC
    LOGICAL :: Found
    INTEGER :: n, nbc
    CHARACTER(*), PARAMETER :: Caller="IntegralProjector"

    nbc = Model % NumberOfBCs
    
    IF(IsBodyForce) THEN    
      BC => Model % BodyForces(BCInd-nbc) % Values
      IF(.NOT. ASSOCIATED(BC)) CALL Warn(Caller,'Why not body force associated!?')
    ELSE
      BC => Model % BCs(BCInd) % Values
    END IF
    NULLIFY(Projector)
    
    IF( .NOT. ListGetLogical( BC,'Integral BC', Found ) ) RETURN
    
    IF(IsBodyForce) THEN
      CALL Info(Caller,'Creating integral constraint matrix for body force: '//I2S(BCind-nbc),Level=6)
    ELSE
      CALL Info(Caller,'Creating integral constraint matrix for boundary: '//I2S(BCind),Level=6)
    END IF
      
    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_INTEGRAL
    
    CALL CreateIntegralProjector()
    
    CALL List_toCRSMatrix(Projector)
    area = SUM( Projector % Values )
    n = SIZE( Projector % Values ) 
    
    WRITE( Message,'(A,ES12.4)') 'Total area of boundary integral:',area  
    CALL Info(Caller, Message, Level=6 )

    CALL SetInvPermIndex()
    
    IF( InfoActive(20) ) THEN
       WRITE(Message,'(A,ES12.3)') 'Sum of constraint matrix entries: ',SUM(Projector%Values)
       CALL Info(Caller,Message)
       CALL Info(Caller,'Constraint matrix cols min: '//I2S(MINVAL(Projector%Cols)))
       CALL Info(Caller,'Constraint matrix cols max: '//I2S(MAXVAL(Projector%Cols)))
       CALL Info(Caller,'Constraint matrix rows min: '//I2S(MINVAL(Projector%Rows)))
       CALL Info(Caller,'Constraint matrix rows max: '//I2S(MINVAL(Projector%Rows)))
     END IF
            
  CONTAINS
    
    SUBROUTINE CreateIntegralProjector()
    
      INTEGER :: i,j,n,t,p,t1,t2
      REAL(KIND=dp) :: u,v,w,weight,x,detJ,val
      REAL(KIND=dp), ALLOCATABLE :: Basis(:)
      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t), POINTER :: Element
      INTEGER, POINTER :: Indexes(:)  
      TYPE(GaussIntegrationPoints_t) :: IP
      LOGICAL :: AxisSym, Stat, Visited = .FALSE.

      SAVE Visited, Nodes, Basis

      IF(.NOT. Visited ) THEN
        n = Mesh % MaxElementNodes
        ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )
        Visited = .TRUE.
      END IF

      AxisSym = ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
          CurrentCoordinateSystem() == CylindricSymmetric ) 

      IF(IsBodyForce) THEN
        t1 = 1
        t2 = Mesh % NumberOfBulkElements 
      ELSE
        t1 = Mesh % NumberOfBulkElements + 1
        t2 = (t1-1) + Mesh % NumberOfBoundaryElements
      END IF
        
      
      DO t = t1, t2

        Element => Mesh % Elements( t )

        IF( IsBodyForce ) THEN
          i = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Stat)
          IF(i /= BCind-nbc) CYCLE
        ELSE
          IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BCInd) % Tag ) CYCLE
        END IF
          
        n = Element % TYPE % NumberOfNodes        
        Indexes => Element % NodeIndexes      
        IP = GaussPoints( Element )

        Nodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
        Nodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))

        DO j=1,IP % n
          u = IP % u(j)
          v = IP % v(j)
          w = IP % w(j)

          Stat = ElementInfo(Element, Nodes, u, v, w, detJ, Basis)

          weight = detJ * IP % s(j)
          IF( AxisSym ) THEN
            x = SUM( Basis(1:n) * Nodes % x(1:n) )
            weight = weight * x
          END IF
          
          DO p=1,n
            val = weight * Basis(p)
            CALL List_AddToMatrixElement(Projector % ListMatrix, 1, Indexes(p), val ) 
          END DO
          
        END DO
      END DO

    END SUBROUTINE CreateIntegralProjector    


    ! Let us associate the inverse permutation to some degree of freedom that is unique and not
    ! set by some other BC / BodyForce. This unique index is needed in the future. 
    !------------------------------------------------------------------------------------------
    SUBROUTINE SetInvPermIndex()
    
      INTEGER :: i,j,t,t1,t2,n,maxind
      TYPE(Element_t), POINTER :: Element
      INTEGER, POINTER :: Indexes(:)  
      LOGICAL, ALLOCATABLE :: SomeOtherBC(:)
      
      IF(.NOT. ASSOCIATED( Projector % InvPerm ) ) THEN
        ALLOCATE( Projector % InvPerm(1) ) 
        Projector % InvPerm = 0
      END IF

      n = Mesh % NumberOfNodes
      ALLOCATE( SomeOtherBC(n) )
      SomeOtherBC = .FALSE.
      maxind = 0
      
      IF(IsBodyForce) THEN
        t1 = 1
        t2 = Mesh % NumberOfBulkElements 
      ELSE
        t1 = Mesh % NumberOfBulkElements + 1
        t2 = (t1-1) + Mesh % NumberOfBoundaryElements
      END IF
        
      DO t = t1, t2
        Element => Mesh % Elements(t)
        IF( IsBodyForce ) THEN
          i = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found)
          IF(i == BCind-nbc) CYCLE
        ELSE          
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(BCInd) % Tag ) CYCLE
        END IF
        Indexes => Element % NodeIndexes      
        SomeOtherBC(Indexes) = .TRUE.
      END DO

      DO t = t1, t2
        Element => Mesh % Elements(t)

        IF( IsBodyForce ) THEN
          i = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force',Found)
          IF(i /= BCind-nbc) CYCLE
        ELSE
          IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BCInd) % Tag ) CYCLE
        END IF
          
        Indexes => Element % NodeIndexes      
        n = Element % TYPE % NumberOfNodes        
        DO i=1,n
          j = Indexes(i)
          IF( SomeOtherBC(j) ) CYCLE
          maxind = MAX(maxind,j)
        END DO
      END DO
      
      IF( maxind == 0 ) THEN
        CALL Fatal(Caller,'Could not determine maximum unset index!')
      ELSE
        CALL Info(Caller,'Setting the representative node to: '//I2S(maxind),Level=8)
        Projector % InvPerm(1) = maxind
      END IF        
    END SUBROUTINE SetInvPermIndex
    
  END FUNCTION IntegralProjector


  


  ! Routine for increasing element order by adding an additional node an each edge.
  ! Basically the same order of elements could be created by p-elements but this provides
  ! alternative solution when nodal finite element are preferred. Often the mesh may be
  ! made quadratic with the preprocessors but this enables also the use of mesh extrusion
  ! and mesh multiplication which cannot be used with higher order nodal elements.
  !--------------------------------------------------------------------------------------
  SUBROUTINE IncreaseElementOrder( Model, Mesh )
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), TARGET :: Mesh
    TYPE(Element_t), POINTER :: Element, Edge
    INTEGER :: n0,n1,m1,m2,i,i1,i2,t,ElemType, NewType, Tinds(4)
    INTEGER, POINTER  :: NewIndexes(:)
    REAL(KIND=dp), POINTER :: x(:), y(:), z(:), xtmp(:)
    
    CALL Info('IncreaseElementOrder','Increasing element order from linear to quadratic!')
    
    IF ( .NOT.ASSOCIATED( Mesh % Edges ) ) THEN
      CALL FindMeshEdges( Mesh )
    END IF
      
    n0 = Mesh % NumberOfNodes
    n1 = Mesh % NumberOfEdges

    CALL Info('IncreaseElementOrder','Adding node to each edge: '//I2S(n1),Level=8)
    
    ! Increase size of coordinate vectors
    ALLOCATE(xtmp(n0))
    xtmp = Mesh % Nodes % x
    DEALLOCATE( Mesh % Nodes % x)
    ALLOCATE( Mesh % Nodes % x(n0+n1))
    x => Mesh % Nodes % x
    x(1:n0) = xtmp; x(n0+1:n0+n1) = 0.0_dp

    xtmp = Mesh % Nodes % y
    DEALLOCATE( Mesh % Nodes % y)
    ALLOCATE( Mesh % Nodes % y(n0+n1))
    y => Mesh % Nodes % y
    y(1:n0) = xtmp; y(n0+1:n0+n1) = 0.0_dp

    xtmp = Mesh % Nodes % z
    DEALLOCATE( Mesh % Nodes % z)
    ALLOCATE( Mesh % Nodes % z(n0+n1))
    z => Mesh % Nodes % z
    z(1:n0) = xtmp; z(n0+1:n0+n1) = 0.0_dp
    DEALLOCATE(xtmp)

    ! Locate new nodes at the center of edges
    DO i=1,Mesh % NumberOfEdges
      Edge => Mesh % Edges(i)
      i1 = Edge % NodeIndexes(1)
      i2 = Edge % NodeIndexes(2)
      x(n0+i) = 0.5_dp*(x(i1)+x(i2))
      y(n0+i) = 0.5_dp*(y(i1)+y(i2))
      z(n0+i) = 0.5_dp*(z(i1)+z(i2))
    END DO

    ! Add the new nodes to the linear elements and
    ! change the element type to reflect the increase in number of nodes.
    DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      ElemType = Element % TYPE % ElementCode
      IF( ElemType == 101) CYCLE
      
      SELECT CASE( ElemType )
      CASE( 101 )
        CYCLE
      CASE( 202, 303, 404, 504, 605, 706, 808 )              
        m1 = Element % TYPE % NumberOfNodes
        m2 = Element % TYPE % NumberOfEdges
        NewType = ElemType + m2
        ALLOCATE( NewIndexes(m1+m2) )
        NewIndexes(1:m1) = Element % NodeIndexes(1:m1)
        NewIndexes(m1+1:m1+m2) = n0 + Element % EdgeIndexes(1:m2)      
        
        IF( ElemType == 808 ) THEN
          ! This is somewhat annoying that the edges and nodes cannot be consistent...
          Tinds(1:4) = NewIndexes(17:20)
          NewIndexes(17:20) = NewIndexes(13:16)
          NewIndexes(13:16) = Tinds(1:4)
        END IF

        DEALLOCATE( Element % NodeIndexes )
        Element % NodeIndexes => NewIndexes
        NULLIFY(NewIndexes)
        Element % TYPE => GetElementType( NewType ) 
      CASE DEFAULT
        CALL Fatal('IncreaseElementOrder','Cannot increase element order for: '//I2S(ElemType))
      END SELECT

    END DO

    ! Parallel info is needed to renumber the nodes in parallel.
    CALL IncreaseParallelInfoOrder()
    
    Mesh % NumberOfNodes = n0 + n1

    CALL ReleaseMeshEdgeTables( Mesh )
    CALL ReleaseMeshFaceTables( Mesh )     

    CALL Info('IncreaseElementOrder','Elements increased to 2nd order serendipity elements')
    
    
  CONTAINS

    
    SUBROUTINE IncreaseParallelInfoOrder()
      TYPE( ParallelInfo_t), POINTER :: ParInfo
      INTEGER, POINTER :: globaldofs(:)
      LOGICAL, POINTER :: ginterface(:)
      TYPE(NeighbourList_t), POINTER  :: NeighbourList(:)
      INTEGER :: globaln0 
      
      IF(ParEnv % PEs == 1 .OR. Mesh % SingleMesh ) RETURN

      ParInfo => Mesh % ParallelInfo
      
      ginterface => ParInfo % Ginterface
      NULLIFY( ParInfo % Ginterface)
      ALLOCATE( ParInfo % Ginterface(n0+n1))
      ParInfo % Ginterface(1:n0) = ginterface(1:n0)
      DEALLOCATE(ginterface)

      globaldofs => ParInfo % Globaldofs
      NULLIFY( ParInfo % Globaldofs)
      ALLOCATE( ParInfo % Globaldofs(n0+n1))
      ParInfo % Globaldofs(1:n0) = globaldofs(1:n0)

      globaln0 = MAXVAL( globaldofs(1:n0) )
      globaln0 = ParallelReduction(globaln0,2)

      DEALLOCATE(globaldofs)
      DO i=1,n1
        ParInfo % Globaldofs(n0+i) = globaln0 + Mesh % Edges(i) % GelementIndex 
      END DO
      
      neighbourList => ParInfo % NeighbourList
      NULLIFY( ParInfo % NeighbourList )
      ALLOCATE( ParInfo % NeighbourList(n0+n1))
      DO i=1, n0
        ParInfo % NeighbourList(i) % Neighbours => NeighbourList(i) % Neighbours
        NULLIFY( NeighbourList(i) % Neighbours )
      END DO
      DEALLOCATE( NeighbourList )

      DO i=1,n1
        ParInfo % NeighbourList(n0+i) % Neighbours => ParInfo % EdgeNeighbourList(i) % Neighbours
        NULLIFY( ParInfo % EdgeNeighbourList(i) % Neighbours )
      END DO

    END SUBROUTINE IncreaseParallelInfoOrder
          
  END SUBROUTINE IncreaseElementOrder


  


!------------------------------------------------------------------------------
!> Generate element edge (faces in 3D) tables for given mesh.
!> Currently only for triangles and tetras. If mesh already
!> has edges do nothing.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges( Mesh, FindEdges, FindFaces )
!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     LOGICAL, OPTIONAL :: FindEdges, FindFaces

     LOGICAL :: FindEdges3D, FindFaces3d
     INTEGER :: MeshDim, SpaceDim, MaxElemDim 

     IF(PRESENT(FindEdges)) THEN
       FindEdges3D = FindEdges
     ELSE
       FindEdges3D = .TRUE.
     END IF

     IF(PRESENT(FindFaces)) THEN
       FindFaces3D = FindFaces
     ELSE
       FindFaces3D = .TRUE.
     END IF

!------------------------------------------------------------------------------

     SpaceDim = CoordinateSystemDimension()
     MeshDim = Mesh % MeshDim

     IF( MeshDim == 0 ) THEN
       CALL Fatal('FindMeshEdges','Mesh dimension is zero!')
     END IF
     IF( SpaceDim > MeshDim ) THEN
       CALL Warn('FindMeshEdges','Mesh dimension and space dimension differ: '&
           // I2S(MeshDim)//' vs. '//I2S(SpaceDim))
     END IF

     MaxElemDim = EnsureElemDim( MeshDim ) 
     IF( MaxElemDim < MeshDim ) THEN
       CALL Warn('FindMeshEdges','Element dimension smaller than mesh dimension: '//&
           I2S(MaxElemDim)//' vs '//I2S(MeshDim))
     END IF


     SELECT CASE( MaxElemDim )

     CASE(1)
       IF ( .NOT.ASSOCIATED( Mesh % Edges ) ) THEN
         CALL Info('FindMeshEdges','Determining edges in 1D mesh',Level=10)
         CALL FindMeshEdges2D( Mesh )
       END IF

     CASE(2)
       IF ( .NOT.ASSOCIATED( Mesh % Edges ) ) THEN
         CALL Info('FindMeshEdges','Determining edges in 2D mesh',Level=10)
         CALL FindMeshEdges2D( Mesh )
       END IF

     CASE(3)
       IF ( .NOT.ASSOCIATED(Mesh % Faces) .AND. FindFaces3D ) THEN
         CALL Info('FindMeshEdges','Determining faces in 3D mesh',Level=10)
         CALL FindMeshFaces3D( Mesh )
       END IF
       IF(FindEdges3D) THEN
         IF ( .NOT.ASSOCIATED( Mesh % Edges) ) THEN
           CALL Info('FindMeshEdges','Determining edges in 3D mesh',Level=10)
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
          IF(FaceInd(j)<=0) CYCLE

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
  SUBROUTINE FindMeshEdges2D( Mesh, BulkMask )
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: BulkMask(:)
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

    LOGICAL :: Found,Masked, LG
    INTEGER :: i,j,k,n,NofEdges,Edge,Swap,Node1,Node2,istat,Degree,maxedges,allocstat
!------------------------------------------------------------------------------
!
!   Initialize:
!   -----------

    CALL Info('FindMeshEdges2D','Finding mesh edges in 2D mesh',Level=12)
    
    Masked = PRESENT(BulkMask)
    
    DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       IF(.NOT.ASSOCIATED(Element)) CYCLE
       IF(Element % Type % ElementCode < 200) CYCLE

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)

           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
               j=Element % Boundaryinfo % Right % ElementIndex
           END IF

           IF(j==-1) CYCLE
         END IF
         IF ( .NOT.BulkMask(j)) CYCLE
       END IF

       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector( Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    CALL Info('FindMeshEdges2D','Creating hash table of size '&
        //I2S(Mesh % NumberOfNodes)//' for node-to-node connectivity',Level=20)
    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
      NULLIFY( HashTable(i) % Head )
    END DO
    CALL Info('FindMeshEdges2D','Hash table allocated',Level=25)
     
!------------------------------------------------------------------------------

    Edges => NULL()
    NofEdges = 0
1   DO i=1,Mesh % NumberOfBulkELements+Mesh % NumberOfBoundaryElements

       Element => Mesh % Elements(i)

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
               j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)
           
           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
                 j=Element % Boundaryinfo % Right % ElementIndex
           END IF
           
           IF(j==-1) CYCLE
         END IF
         
         IF(.NOT. BulkMask(j)) CYCLE
       END IF

       SELECT CASE( Element % TYPE % ElementCode / 100 )
       CASE(1) 
         CYCLE
       CASE(2)
         n = 1
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
         IF(n==1) THEN
           Node2 = Element % NodeIndexes(2)
         ELSE IF ( k<n ) THEN
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

         IF(.NOT. ASSOCIATED( Edges ) ) THEN
           ! Edge has already been numbered
           IF(Found ) CYCLE

           ! This is visited only the first round when Edges have not been allocated.           
           NofEdges = NofEdges + 1
           Edge = NofEdges
           
           ! Update the hash table:
           !----------------------
           ALLOCATE( HashPtr, STAT=allocstat )
           IF( allocstat /= 0 ) THEN
             CALL Fatal('FindMeshEdges2D','Allocation error for HashPtr allocation')
           END IF           
           HashPtr % Edge = Edge
           HashPtr % Node = Node2
           HashPtr % Next => HashTable(Node1) % Head
           HashTable(Node1) % Head => HashPtr
         
         ELSE 
           IF(.NOT. Found ) THEN
             CALL Fatal('FindMeshEdges2D','We should find the edge in the hash table!')
           END IF
           IF( Edge > SIZE( Edges ) ) THEN
             CALL Fatal('FindMeshEdges2D','Number of edges larger than expected!')
           END IF
                      
           IF(.NOT. ASSOCIATED(Edges(Edge) % TYPE ) ) THEN
             Degree = MAX( Element % TYPE % BasisFunctionDegree, 1)

             Edges(Edge) % ElementIndex = Edge
             CALL AllocateVector( Edges(Edge) % NodeIndexes, Degree+1)
             ALLOCATE( Edges(Edge) % BoundaryInfo, STAT=allocstat )
             IF( allocstat /= 0 ) THEN
               CALL Fatal('FindMeshEdges2D','Allocation error for BoyndaryInfo allocation')
             END IF
             Edges(Edge) % TYPE => GetElementType( 201+Degree, .FALSE. )

             Edges(Edge) % NodeIndexes(1) = Element % NodeIndexes(k)
             IF( n==1 ) THEN
               Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(2)
             ELSE IF ( k < n ) THEN
               Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(k+1)
             ELSE
               Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(1)
             END IF

             DO j=2,Degree
               Edges(Edge) % NodeIndexes(j+1) = Element % NodeIndexes(k+n+j-2)
             END DO
             Edges(Edge) % PartIndex = Element % PartIndex
             
             ! Create P element definitions if needed
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
               CALL AllocatePDefinitions(Edges(Edge))
               Edges(Edge) % PDefs % P = 0
             ELSE
               NULLIFY( Edges(Edge) % PDefs )
             END IF

             Edges(Edge) % NDofs = 0
             IF (Element % NDOFs /= 0 ) Edges(Edge) % NDOFs = &
                 Element % NDOFs / Element % TYPE % NumberOfNodes * &
                 Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             NULLIFY( Edges(Edge) % EdgeIndexes )
             NULLIFY( Edges(Edge) % FaceIndexes )
             
             Edges(Edge) % BoundaryInfo % Left  => NULL()
             Edges(Edge) % BoundaryInfo % Right => NULL()
           END IF

           ! These structures need to be updated to both new and old edge.
           Element % EdgeIndexes(k) = Edge
           IF (i <= Mesh % NumberofBulkElements) THEN
             IF(ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
               Edges(Edge) % BoundaryInfo % Right => Element
             ELSE
               Edges(Edge) % BoundaryInfo % Left => Element
             END IF
           END IF
           
         END IF
       END DO
     END DO

     IF(.NOT. ASSOCIATED( Edges ) ) THEN
       CALL Info('FindMeshEdges2D','Allocating edge table of size: '//I2S(NofEdges),Level=12)
       CALL AllocateVector( Mesh % Edges, NofEdges ) 
       Edges => Mesh % Edges
       GOTO 1
     END IF
         
    Mesh % NumberOfEdges = NofEdges
    CALL Info('FindMeshEdges2D','Number of edges found: '//I2S(NofEdges),Level=10)

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

    CALL Info('FindMeshEdges2D','All done',Level=20)

!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh faces.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshFaces3D( Mesh, BulkMask)
    USE PElementMaps, ONLY : GetElementFaceMap

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: BulkMask(:)
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

    LOGICAL :: Found,Masked,LG
    INTEGER :: n1,n2,n3,n4
    INTEGER :: i,j,k,n,NofFaces,Face,Swap,Node1,Node2,Node3,istat,Degree,facenodes
     
    TYPE(Element_t), POINTER :: Element, Faces(:)

    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET  :: TetraFaceMap(4,6), BrickFaceMap(6,9), &
         WedgeFaceMap(5,8), PyramidFaceMap(5,8), TriFaceMap(1,3), QuadFaceMap(1,4)
    
    INTEGER :: nf(4)
!------------------------------------------------------------------------------
    
    CALL Info('FindMeshFaces3D','Finding mesh faces in 3D mesh',Level=12)

    Masked = PRESENT(BulkMask)

    TriFaceMap(1,:)  = [1,2,3]
    QuadFaceMap(1,:) = [1,2,3,4]

    TetraFaceMap(1,:) = [ 1, 2, 3, 5, 6, 7 ]
    TetraFaceMap(2,:) = [ 1, 2, 4, 5, 9, 8 ]
    TetraFaceMap(3,:) = [ 2, 3, 4, 6, 10, 9 ]
    TetraFaceMap(4,:) = [ 3, 1, 4, 7, 8,10 ]

    WedgeFaceMap(1,:) = [ 1, 2, 3, 7, 8, 9, -1, -1 ]
    WedgeFaceMap(2,:) = [ 4, 5, 6, 10, 11, 12, -1, -1 ]
    WedgeFaceMap(3,:) = [ 1, 2, 5, 4, 7, 14, 10, 13 ]
    WedgeFaceMap(4,:) = [ 3, 2, 5, 6, 8, 14, 11, 15 ]
    WedgeFaceMap(5,:) = [ 3, 1, 4, 6, 9, 13, 12, 15 ]

    PyramidFaceMap(1,:) = [ 1, 2, 3, 4,  6,  7,  8,  9 ]
    PyramidFaceMap(2,:) = [ 1, 2, 5, 6, 11, 10, -1, -1 ]
    PyramidFaceMap(3,:) = [ 2, 3, 5, 7, 12, 11, -1, -1 ]
    PyramidFaceMap(4,:) = [ 3, 4, 5, 8, 13, 12, -1, -1 ]
    PyramidFaceMap(5,:) = [ 4, 1, 5, 9, 10, 13, -1, -1 ]

    BrickFaceMap(1,:) = [ 1, 2, 3, 4,  9, 10, 11, 12, 25 ]
    BrickFaceMap(2,:) = [ 5, 6, 7, 8, 17, 18, 19, 20, 26 ]
    BrickFaceMap(3,:) = [ 1, 2, 6, 5,  9, 14, 17, 13, 21 ]
    BrickFaceMap(4,:) = [ 2, 3, 7, 6, 10, 15, 18, 14, 22 ]
    BrickFaceMap(5,:) = [ 3, 4, 8, 7, 11, 16, 19, 15, 23 ]
    BrickFaceMap(6,:) = [ 4, 1, 5, 8, 12, 13, 20, 16, 24 ]

!
!   Initialize:
!   -----------   
    DO i=1,SIZE(Mesh % Elements)
       Element => Mesh % Elements(i)

       IF(.NOT.ASSOCIATED(Element % Type)) CYCLE
       IF(Element % Type % ElementCode<300 ) CYCLE

       IF(Masked) THEN
         j = i
         IF(i>Mesh % NumberOfBulkElements) THEN
           j = -1
           IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

           LG=.FALSE.
           IF(j>0) LG=BulkMask(j)

           IF(.NOT. LG) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
               j=Element % Boundaryinfo % Right % ElementIndex
           END IF

           IF(j==-1) CYCLE
         END IF

         IF(.NOT. BulkMask(j)) CYCLE
       END IF

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
    Faces => NULL()

1   DO i=1,SIZE(Mesh % Elements)

      Element => Mesh % Elements(i)
      IF(.NOT.ASSOCIATED(Element % Type)) CYCLE
      IF(Element % Type % ElementCode < 300 ) Cycle

      IF(Masked) THEN
        j = i
        IF(i>Mesh % NumberOfBulkElements) THEN
          j = -1
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) &
              j=Element % Boundaryinfo % Left % ElementIndex

          LG=.FALSE.
          IF(j>0) LG=BulkMask(j)

          IF(.NOT. LG) THEN
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) &
                j=Element % Boundaryinfo % Right % ElementIndex
          END IF

          IF(j==-1) CYCLE
        END IF
        IF(.NOT. BulkMask(j)) CYCLE
      END IF

      ! For P elements mappings are different
      IF ( ASSOCIATED(Element % PDefs) ) THEN
        CALL GetElementFaceMap(Element, FaceMap)
        n = Element % TYPE % NumberOfFaces
      ELSE
        SELECT CASE( Element % TYPE % ElementCode / 100 )
        CASE(3)
          n = 1
          FaceMap => TriFaceMap
        CASE(4)
          n = 1
          FaceMap => QuadFaceMap
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
          CALL Fatal('FindMeshFaces','Element type '&
              //I2S(Element % Type % ElementCode)//' not implemented!')
        END SELECT
      END IF
 
!      Loop over every face of every element:
!      --------------------------------------
      DO k=1,n
                    
        SELECT CASE( Element % TYPE % ElementCode / 100 )
          
        CASE(3)
          ! Triangle:
          !=======
          facenodes = 3

        CASE(4)
          ! Quad:
          !=======
          facenodes = 4

        CASE(5)
          ! Tetras:
          !=======
          facenodes = 3

        CASE(6)
          ! Pyramids:
          !=========
          IF ( k == 1 ) THEN
            facenodes = 4
          ELSE
            facenodes = 3
          END IF
          
        CASE(7)
          ! Wedges:
          !=======
          IF ( k <= 2 ) THEN
            facenodes = 3
          ELSE
            facenodes = 4
          END IF
                
        CASE(8)
          ! Bricks:
          !=======
          facenodes = 4
          
        CASE DEFAULT
          WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
          CALL Fatal('FindMeshFaces',Message)
        END SELECT

        nf(1:facenodes) = Element % NodeIndexes(FaceMap(k,1:facenodes))
        CALL sort( facenodes, nf )
        
!         We use MIN(Node1,Node2,Node3) as the hash table key:
!         ---------------------------------------------------
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
        
!         Existing face, update structures:
!         ----------------------------------

        IF( .NOT. ASSOCIATED( Faces ) ) THEN
          IF(Found ) CYCLE

          ! Update the hash table:
          !----------------------
          NofFaces = NofFaces + 1
          Face = NofFaces
          ALLOCATE( HashPtr )
          HashPtr % Face = Face
          HashPtr % Node1 = Node2
          HashPtr % Node2 = Node3
          HashPtr % Next => HashTable(Node1) % Head
          HashTable(Node1) % Head => HashPtr
        ELSE
          IF(.NOT. Found ) THEN
            CALL Fatal('FindMeshFaces3D','We should find the edge in the hash table!')
          END IF
          IF( Face > SIZE( Faces ) ) THEN
            CALL Fatal('FindMeshFaces3D','Number of faces larger than expected!')
          END IF
          
          IF(.NOT. ASSOCIATED( Faces(Face) % TYPE ) ) THEN
            ! Face not yet there, create:
            !---------------------------
            Degree = Element % TYPE % BasisFunctionDegree
            Faces(Face) % ElementIndex = Face
            
            SELECT CASE( Element % TYPE % ElementCode / 100 )

            CASE(1,2)
              CYCLE

            CASE(3)
              ! linear tri
              !-----------
              SELECT CASE( Degree ) 
              CASE(1)
                n1 = 3
              CASE DEFAULT
              END SELECT
              Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              
            CASE(4)
              ! linear quad
              !-----------
              SELECT CASE( Degree ) 
              CASE(1)
                n1 = 4
              CASE DEFAULT
              END SELECT              
              Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
              
            CASE(5)
              ! for tetras:
              !-----------
              SELECT CASE( Degree ) 
              CASE(1)
                n1 = 3
              CASE(2)
                n1 = 6
              CASE(3)
                n1 = 10
              END SELECT
              
              Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              
            CASE(6)              
               ! Pyramids ( 605 and 613 supported )
               !-------------------------------
              IF ( k == 1 ) THEN
                n1 = Degree * 4
                Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
              ELSE
                n1 = Degree * 3
                Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              END IF
              
            CASE(7)
               ! for wedges, 706 and 715 supported:
               !-------------------------------
              IF ( k <= 2 ) THEN
                n1 = Degree * 3
                Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
              ELSE
                n1 = Degree * 4
                Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
              END IF
              
            CASE(8)
               ! for bricks:
               !-----------
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
              CALL Fatal('FindMeshFaces','Element type '&
                  //I2S(Element % TYPE % ElementCode)//' not implemented!')
              
            END SELECT
            
             ! Allocate p structures for p elements
            IF ( ASSOCIATED( Element % PDefs ) ) THEN
              CALL AllocatePDefinitions(Faces(Face))
              Faces(Face) % PDefs % P = 0
            ELSE
              NULLIFY( Faces(Face) % PDefs )
            END IF
            
            Faces(Face) % NDOFs  = 0
            IF (Element % NDOFs /= 0) Faces(Face) % NDOFs = &
                Element % NDOFs / Element % TYPE % NumberOfNodes * &
                Faces(Face) % TYPE % NumberOfNodes
            Faces(Face) % BDOFs  = 0
            Faces(Face) % DGDOFs = 0
            Faces(Face) % EdgeIndexes => NULL()
            Faces(Face) % FaceIndexes => NULL()
            
            CALL AllocateVector( Faces(Face) % NodeIndexes,n1 )
            DO n2=1,n1
              Faces(Face) % NodeIndexes(n2) = &
                  Element % NodeIndexes(FaceMap(k,n2)) 
            END DO

            Faces(Face) % PartIndex = Element % PartIndex

            ALLOCATE( Faces(Face) % BoundaryInfo )
            Faces(Face) % BoundaryInfo % Left  => NULL()
            Faces(Face) % BoundaryInfo % Right => NULL()
          END IF

          Element % FaceIndexes(k) = Face            
          IF(i<=Mesh % NumberOfBulkElements) THEN
            IF( ASSOCIATED(Faces(Face) % BoundaryInfo % Left) ) THEN
              Faces(Face) % BoundaryInfo % Right => Element
            ELSE
              Faces(Face) % BoundaryInfo % Left => Element
            END IF
          END IF
          
        END IF
      END DO
    END DO

    IF(.NOT. ASSOCIATED( Faces ) ) THEN
      CALL Info('FindMeshFaces3D','Allocating face table of size: '&
          //I2S(NofFaces),Level=25)
      CALL AllocateVector( Mesh % Faces, NofFaces, 'FindMeshFaces3D' )
      Faces => Mesh % Faces
      GOTO 1
    END IF
        
    Mesh % NumberOfFaces = NofFaces
    CALL Info('FindMeshFaces3D','Number of faces found: '//I2S(NofFaces),Level=10)

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

    CALL Info('FindMeshFaces3D','All done',Level=20)
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshFaces3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh edges.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges3D( Mesh )
    USE PElementMaps, ONLY : GetElementEdgeMap, GetElementFaceEdgeMap

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
    INTEGER :: n1,n2, n_e, maxedges
    INTEGER :: i,j,k,n,NofEdges,Edge,Node1,Node2,istat,Degree,ii,jj
     
    TYPE(Element_t), POINTER :: Element, Edges(:), Face

    INTEGER, POINTER :: EdgeMap(:,:), FaceEdgeMap(:,:)
    INTEGER, TARGET  :: TetraEdgeMap(6,3), BrickEdgeMap(12,3), TetraFaceMap(4,6), &
      WedgeEdgeMap(9,3), PyramidEdgeMap(8,3), TetraFaceEdgeMap(4,3), &
      BrickFaceEdgeMap(8,4), WedgeFaceEdgeMap(6,4), PyramidFaceEdgeMap(5,4), &
         QuadEdgeMap(4,3), TriEdgeMap(3,3), TriFaceMap(1,3), QuadFaceMap(1,4), LineEdgeMap(1,2)
!------------------------------------------------------------------------------
    
    CALL Info('FindMeshEdges3D','Finding mesh edges in 3D mesh',Level=12)

    LineEdgeMap(1,:) = [1,2]

    TriEdgeMap(1,:) = [1,2,4]
    TriEdgeMap(2,:) = [2,3,5]
    TriEdgeMap(3,:) = [3,1,6]

    TriFaceMap(1,:) = [1,2,3]

    QuadEdgeMap(1,:) = [1,2,5]
    QuadEdgeMap(2,:) = [2,3,6]
    QuadEdgeMap(3,:) = [3,4,7]
    QuadEdgeMap(4,:) = [4,1,8]

    QuadFaceMap(1,:) = [1,2,3,4]

    TetraFaceMap(1,:) = [ 1, 2, 3, 5, 6, 7 ]
    TetraFaceMap(2,:) = [ 1, 2, 4, 5, 9, 8 ]
    TetraFaceMap(3,:) = [ 2, 3, 4, 6,10, 9 ]
    TetraFaceMap(4,:) = [ 3, 1, 4, 7, 8,10 ]

    TetraFaceEdgeMap(1,:) = [ 1,2,3 ]
    TetraFaceEdgeMap(2,:) = [ 1,5,4 ]
    TetraFaceEdgeMap(3,:) = [ 2,6,5 ]
    TetraFaceEdgeMap(4,:) = [ 3,4,6 ]

    TetraEdgeMap(1,:) = [ 1,2,5 ]
    TetraEdgeMap(2,:) = [ 2,3,6 ]
    TetraEdgeMap(3,:) = [ 3,1,7 ]
    TetraEdgeMap(4,:) = [ 1,4,8 ]
    TetraEdgeMap(5,:) = [ 2,4,9 ]
    TetraEdgeMap(6,:) = [ 3,4,10 ]

    PyramidEdgeMap(1,:) = [ 1,2,1 ]
    PyramidEdgeMap(2,:) = [ 2,3,1 ]
    PyramidEdgeMap(3,:) = [ 3,4,1 ]
    PyramidEdgeMap(4,:) = [ 4,1,1 ]
    PyramidEdgeMap(5,:) = [ 1,5,1 ]
    PyramidEdgeMap(6,:) = [ 2,5,1 ]
    PyramidEdgeMap(7,:) = [ 3,5,1 ]
    PyramidEdgeMap(8,:) = [ 4,5,1 ]

    PyramidFaceEdgeMap(1,:) = [ 1,2,3,4 ]
    PyramidFaceEdgeMap(2,:) = [ 1,6,5,0 ]
    PyramidFaceEdgeMap(3,:) = [ 2,7,6,0 ]
    PyramidFaceEdgeMap(4,:) = [ 3,8,7,0 ]
    PyramidFaceEdgeMap(5,:) = [ 4,5,8,0 ]

    WedgeEdgeMap(1,:) = [ 1, 2, 1 ]
    WedgeEdgeMap(2,:) = [ 2, 3, 1 ]
    WedgeEdgeMap(3,:) = [ 1, 3, 1 ]
    WedgeEdgeMap(4,:) = [ 4, 5, 1 ]
    WedgeEdgeMap(5,:) = [ 5, 6, 1 ]
    WedgeEdgeMap(6,:) = [ 6, 4, 1 ]
    WedgeEdgeMap(7,:) = [ 1, 4, 1 ]
    WedgeEdgeMap(8,:) = [ 2, 5, 1 ]
    WedgeEdgeMap(9,:) = [ 3, 6, 1 ]

    WedgeFaceEdgeMap(1,:) = [ 1,2,3,0 ]
    WedgeFaceEdgeMap(2,:) = [ 4,5,6,0 ]
    WedgeFaceEdgeMap(3,:) = [ 1,8,4,7 ]
    WedgeFaceEdgeMap(4,:) = [ 2,9,5,8 ]
    WedgeFaceEdgeMap(5,:) = [ 3,7,6,9 ]

    BrickEdgeMap(1,:) = [ 1, 2,  9 ]
    BrickEdgeMap(2,:) = [ 2, 3,  10 ]
    BrickEdgeMap(3,:) = [ 4, 3,  11 ]
    BrickEdgeMap(4,:) = [ 1, 4,  12 ]
    BrickEdgeMap(5,:) = [ 5, 6,  13 ]
    BrickEdgeMap(6,:) = [ 6, 7,  14 ]
    BrickEdgeMap(7,:) = [ 8, 7,  15 ]
    BrickEdgeMap(8,:) = [ 5, 8,  16 ]
    BrickEdgeMap(9,:) = [ 1, 5,  17 ]
    BrickEdgeMap(10,:) = [ 2, 6, 18 ]
    BrickEdgeMap(11,:) = [ 3, 7, 19 ]
    BrickEdgeMap(12,:) = [ 4, 8, 20 ]

    BrickFaceEdgeMap(1,:) = [ 1,2,3,4   ]
    BrickFaceEdgeMap(2,:) = [ 5,6,7,8   ]    
    BrickFaceEdgeMap(3,:) = [ 1,10,5,9  ]
    BrickFaceEdgeMap(4,:) = [ 2,11,6,10 ]
    BrickFaceEdgeMap(5,:) = [ 3,12,7,11 ]
    BrickFaceEdgeMap(6,:) = [ 4,9,8,12  ]

!
!   Initialize:
    !   -----------
    n_e = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

    DO i=1,n_e
       Element => Mesh % Elements(i)
       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector(Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    CALL Info('FindMeshEdges3D','Hash table allocated',Level=25)

    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

    !   Loop over elements:
    !   -------------------
    NofEdges = 0
    Edges => NULL()
    
1   DO i=1,n_e
      Element => Mesh % Elements(i)
      
      ! For P elements mappings are different
      IF ( ASSOCIATED(Element % PDefs) ) THEN
        CALL GetElementEdgeMap( Element, EdgeMap )
        IF(Element % Type % ElementCode >= 500) &
          CALL GetElementFaceEdgeMap( Element, FaceEdgeMap ) 

        n = Element % TYPE % NumberOfEdges
      ELSE 
        SELECT CASE( Element % TYPE % ElementCode / 100 )
        CASE(1)
          CYCLE
        CASE(2)
          n = 1
          EdgeMap => LineEdgeMap
          FaceEdgeMap => NULL()
        CASE(3)
          n = 3
          EdgeMap => TriEdgeMap
          FaceEdgeMap => NULL()
        CASE(4)
          n = 4
          EdgeMap => QuadEdgeMap
          FaceEdgeMap => NULL()
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
          CALL Fatal('FindMeshEdges3D','Element type '//I2S(Element % TYPE % ElementCode)//' not implemented!') 
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

        ! Look the edge from the hash table:
        !----------------------------------
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
        
        IF(.NOT. ASSOCIATED( Edges ) ) THEN
          IF( Found ) CYCLE

          NofEdges = NofEdges + 1
          Edge = NofEdges
          
          ! Update the hash table:
          !----------------------
          ALLOCATE( HashPtr )
          HashPtr % Edge = Edge
          HashPtr % Node1 = Node2
          HashPtr % Next => HashTable(Node1) % Head
          HashTable(Node1) % Head => HashPtr
        ELSE
          IF(.NOT. Found ) THEN
            CALL Fatal('FindMeshEdges3D','We should find the edge in the hash table!')
          END IF
          IF( Edge > SIZE( Edges ) ) THEN
            CALL Fatal('FindMeshEdges3D','Number of edges larger than expected!')
          END IF

          Edges(Edge) % ElementIndex = Edge
                    
          IF( ASSOCIATED( Edges(Edge) % TYPE ) ) THEN
            IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
              Edges(Edge) % BoundaryInfo % Left  => Element
            ELSE
              Edges(Edge) % BoundaryInfo % Right => Element
            END IF
          ELSE
            Degree = Element % TYPE % BasisFunctionDegree

            ! Edge is always a line segment with deg+1 nodes:
            !-----------------------------------------------
            Edges(Edge) % TYPE => GetElementType( 201 + degree, .FALSE.)

            Edges(Edge) % NDOFs  = 0
            IF (Element % NDOFs /= 0) Edges(Edge) % NDOFs = &
                Element % NDOFs / Element % TYPE % NumberOfNodes * &
                Edges(Edge) % TYPE % NumberOfNodes
            Edges(Edge) % BDOFs  = 0
            Edges(Edge) % DGDOFs = 0
            Edges(Edge) % EdgeIndexes => NULL()
            Edges(Edge) % FaceIndexes => NULL()
            
            CALL AllocateVector( Edges(Edge) % NodeIndexes, degree + 1 )
            DO n2=1,degree+1
              Edges(Edge) % NodeIndexes(n2) = &
                  Element % NodeIndexes(EdgeMap(k,n2))
            END DO
            
            ALLOCATE( Edges(Edge) % BoundaryInfo )
            Edges(Edge) % BoundaryInfo % Left  => NULL()
            Edges(Edge) % BoundaryInfo % Right => NULL()
            
            ! Allocate P element definitions 
            IF ( ASSOCIATED( Element % PDefs ) ) THEN
              CALL AllocatePDefinitions(Edges(Edge))              
              Edges(Edge) % PDefs % P = 0
            ELSE
              NULLIFY( Edges(Edge) % PDefs )
            END IF            
          END IF

          ! Stuff for both existing and new edge
          !--------------------------------------
          Element % EdgeIndexes(k) = Edge
          
          IF ( ASSOCIATED(Mesh % Faces) .AND. ASSOCIATED(FaceEdgeMap) ) THEN
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
        END IF
          
      END DO
    END DO

    IF(.NOT. ASSOCIATED( Edges ) ) THEN  
      CALL Info('FindMeshEdges3D','Allocating edge table of size: '//I2S(NofEdges),Level=20)
      CALL AllocateVector( Mesh % Edges, NofEdges ) 
      Edges => Mesh % Edges
      CALL Info('FindMeshEdges3D','Edge table allocated',Level=25)
      GOTO 1
    END IF

    Mesh % NumberOfEdges = NofEdges
    CALL Info('FindMeshEdges3D','Number of edges found: '//I2S(NofEdges),Level=10)
    
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

    CALL Info('FindMeshEdges3D','All done',Level=20)

CONTAINS 

    SUBROUTINE FixFaceEdges()

      INTEGER :: i,j,k,n,swap,edgeind(4),i1(2),i2(2)

      DO i=1,Mesh % NumberOfFaces
        Face => Mesh % Faces(i)
        IF(.NOT.ASSOCIATED(Face % EdgeIndexes)) CYCLE
        n = Face % TYPE % NumberOfEdges
        Edgeind(1:n) = Face % EdgeIndexes(1:n)
        IF(ANY(EdgeInd(1:n)==0)) CYCLE
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
  !> Mark edges that define the geometry.
  !> We first identify potential face elements at interface and create mapping
  !> from edges to these faces. Then we check whether any face pair is beyond
  !> a critical angle. 
  !------------------------------------------------------------------------------
  SUBROUTINE MarkSharpEdges( Mesh, SharpEdge, phi0 )
    TYPE(Mesh_t), TARGET :: Mesh
    LOGICAL, ALLOCATABLE :: SharpEdge(:)
    REAL(KIND=dp) :: phi0

    INTEGER :: t,i,i1,i2,j,n,Sweep
    REAL(KIND=dp) :: cosphi, cosphi0, Normal1(3), Normal2(3)
    INTEGER, ALLOCATABLE :: EdgeUses(:), EdgeToFaceMap(:,:)
    TYPE(Element_t), POINTER :: Face1, Face2
    TYPE(Nodes_t), SAVE :: Nodes1, Nodes2
    
    IF(.NOT. ASSOCIATED(Mesh % Faces)) THEN
      CALL FindMeshFaces3D( Mesh )
    END IF    
    IF(.NOT. ASSOCIATED(Mesh % Edges)) THEN
      CALL FindMeshEdges3D(Mesh)
    END IF

    cosphi0 = COS(pi*phi0/180.0_dp)
    
    IF(.NOT. ALLOCATED(SharpEdge)) THEN
      ALLOCATE(SharpEdge(Mesh % NumberOfEdges))
    END IF
    SharpEdge = .FALSE.
    
    n = Mesh % NumberOfEdges
    CALL Info('MarkSharpEdges','Total number of edges '//I2S(n),Level=10)
    ALLOCATE(EdgeUses(n))
    EdgeUses = 0

    ! First mark those face elements that are at interface of two different bodies,
    ! or at outer interface. Note: this does not work yeat in parallel!
    DO Sweep=0,1    
      DO t=1,Mesh % NumberOfFaces
        Face1 => Mesh % Faces(t)

        IF(.NOT. ASSOCIATED(Face1 % BoundaryInfo)) CYCLE
        i1 = 0; i2 = 0
        IF(ASSOCIATED(Face1 % BoundaryInfo % Left)) THEN
          i1 = Face1 % BoundaryInfo % Left % BodyId
        END IF
        IF(ASSOCIATED(Face1 % BoundaryInfo % Right)) THEN
          i2 = Face1 % BoundaryInfo % Right % BodyId
        END IF
        IF(i1 == i2) CYCLE
        
        IF(Sweep == 0) THEN
          ! At first round only count the appearances.
          EdgeUses(Face1 % EdgeIndexes) = EdgeUses(Face1 % EdgeIndexes) + 1
        ELSE
          ! At second round create the mapping from edges to interface faces.
          DO i=1,Face1 % Type % NumberOfEdges
            j = Face1 % EdgeIndexes(i)
            EdgeUses(j) = EdgeUses(j) + 1
            EdgeToFaceMap(j,EdgeUses(j)) = t            
          END DO          
        END IF          
      END DO

      IF(Sweep==0) THEN
        n = MAXVAL(EdgeUses)
        CALL Info('MarkSharpEdges','Edge associated at max. '//I2S(n)//' interface faces',Level=6)
        ALLOCATE(EdgeToFaceMap(Mesh % NumberOfEdges,n))
        EdgeUses = 0
        EdgeToFaceMap = 0
      END IF      
    END DO

    ! Now compute the angle between normals related to faces sharing the edge.
    DO t=1,Mesh % NumberOfEdges    
      DO i1=1, EdgeUses(t)
        Face1 => Mesh % Faces(EdgeToFaceMap(t,i1))
        CALL CopyElementNodesFromMesh(Nodes1,Mesh,&
            Face1 % TYPE % NumberOfNodes,Face1 % NodeIndexes)
        Normal1 = NormalVector(Face1,Nodes1)
        DO i2=i1+1, EdgeUses(t)
          Face2 => Mesh % Faces(EdgeToFaceMap(t,i2))
          CALL CopyElementNodesFromMesh(Nodes2,Mesh,&
              Face2 % TYPE % NumberOfNodes,Face2 % NodeIndexes)
          Normal2 = NormalVector(Face2,Nodes2)
          
          ! Compare cosphi rather than phi since we save one trigonometric operation. 
          cosphi = ABS(SUM(Normal1 * Normal2))
          IF(cosphi < cosphi0) SharpEdge(t) = .TRUE.
        END DO
      END DO
    END DO
       
    n = COUNT(SharpEdge)
    CALL Info('MarkSharpEdges','Number of sharp edges is '//I2S(n),Level=5)

    DEALLOCATE(EdgeUses,EdgeToFaceMap)

#if 0
    ! For debugging reasons we may want to save the edges. 
    ! plot3(sharp(
    OPEN( 10, FILE = 'sharp_edge.dat' )    
    DO t=1, Mesh % NumberOfEdges
      IF(.NOT. SharpEdge(t)) CYCLE
      i1 = Mesh % Edges(t) % NodeIndexes(1)
      i2 = Mesh % Edges(t) % NodeIndexes(2)
      WRITE(10,*) t,Mesh % Nodes % x(i1),Mesh % Nodes % y(i1),Mesh % Nodes % z(i1), &
          Mesh % Nodes % x(i2),Mesh % Nodes % y(i2),Mesh % Nodes % z(i2)
    END DO
    CLOSE(10)
#endif
    
  END SUBROUTINE MarkSharpEdges


  SUBROUTINE MarkSharpNodes( Mesh, SharpEdge, SharpNode, phi0 )
    TYPE(Mesh_t), TARGET :: Mesh
    LOGICAL :: SharpEdge(:)
    REAL(KIND=dp) :: phi0
    LOGICAL, ALLOCATABLE :: SharpNode(:)

    INTEGER :: t,i,j,i1,i2,j1,j2,n,Sweep
    REAL(KIND=dp) :: cosphi, cosphi0, Normal1(3), Normal2(3)
    INTEGER, ALLOCATABLE :: NodeUses(:), NodeToEdgeMap(:,:)
    TYPE(Element_t), POINTER :: Edge1, Edge2
    
    IF(.NOT. ASSOCIATED(Mesh % Edges)) THEN
      CALL Fatal('MarkSharpNodes','We should have edges allocated!')
    END IF

    cosphi0 = COS(pi*phi0/180.0_dp)
    
    IF(.NOT. ALLOCATED(SharpNode)) THEN
      ALLOCATE(SharpNode(Mesh % NumberOfNodes))
    END IF
    SharpNode = .FALSE.

    n = Mesh % NumberOfNodes
    CALL Info('MarkSharpNodes','Total number of nodes '//I2S(n),Level=10)
    ALLOCATE(NodeUses(n))
    NodeUses = 0
    
    ! First create a structure from potential corner nodes to all sharp edges.
    DO Sweep=0,1    
      DO t=1,Mesh % NumberOfEdges
        IF(.NOT. SharpEdge(t)) CYCLE
                
        Edge1 => Mesh % Edges(t)        
        IF(Sweep == 0) THEN
          ! At first round only count the appearances.
          NodeUses(Edge1 % NodeIndexes) = NodeUses(Edge1 % NodeIndexes) + 1
        ELSE
          ! At second round create the mapping from edges to interface faces.
          DO i=1,Edge1 % Type % NumberOfNodes
            j = Edge1 % NodeIndexes(i)
            NodeUses(j) = NodeUses(j) + 1
            NodeToEdgeMap(j,NodeUses(j)) = t            
          END DO          
        END IF          
      END DO

      IF(Sweep==0) THEN
        n = MAXVAL(NodeUses)
        CALL Info('MarkSharpNodes','Node associated at max. '//I2S(n)//' sharp edges',Level=6)
        ALLOCATE(NodeToEdgeMap(Mesh % NumberOfNodes,n))
        n = COUNT(NodeUses > 1)
        CALL Info('MarkSharpNodes','Number of sharp candidate nodes is '//I2S(n),Level=6)
        NodeUses = 0
        NodeToEdgeMap = 0
      END IF
    END DO

    ! Now compute the angle between edges related to the potential corner node. 
    DO t=1,Mesh % NumberOfNodes    
      DO i1=1, NodeUses(t)
        Edge1 => Mesh % Edges(NodeToEdgeMap(t,i1))
        j1 = Edge1 % NodeIndexes(1)
        j2 = Edge1 % NodeIndexes(2)        
        Normal1(1) = Mesh % Nodes % x(j1) - Mesh % Nodes % x(j2)
        Normal1(2) = Mesh % Nodes % y(j1) - Mesh % Nodes % y(j2)
        Normal1(3) = Mesh % Nodes % z(j1) - Mesh % Nodes % z(j2)
        Normal1 = Normal1 / SQRT(SUM(Normal1*Normal1))
        
        DO i2=i1+1, NodeUses(t)
          Edge2 => Mesh % Edges(NodeToEdgeMap(t,i2))
          j1 = Edge2 % NodeIndexes(1)
          j2 = Edge2 % NodeIndexes(2)        
          Normal2(1) = Mesh % Nodes % x(j1) - Mesh % Nodes % x(j2)
          Normal2(2) = Mesh % Nodes % y(j1) - Mesh % Nodes % y(j2)
          Normal2(3) = Mesh % Nodes % z(j1) - Mesh % Nodes % z(j2)
          Normal2 = Normal2 / SQRT(SUM(Normal2*Normal2))

          ! Compare cosphi rather than phi since we save one trigonometric operation. 
          cosphi = ABS(SUM(Normal1 * Normal2))
          IF(cosphi < cosphi0) SharpNode(t) = .TRUE.
        END DO
      END DO
    END DO
       
    n = COUNT(SharpNode)
    CALL Info('MarkSharpNodes','Number of sharp nodes is '//I2S(n),Level=5)

    DEALLOCATE(NodeUses,NodeToEdgeMap)

#if 0
    ! For debugging reasons we may want to save the corner nodes. 
    OPEN( 10, FILE = 'sharp_node.dat' )    
    DO t=1, Mesh % NumberOfNodes
      IF(.NOT. SharpNode(t)) CYCLE
      WRITE(10,*) t,Mesh % Nodes % x(t),Mesh % Nodes % y(t),Mesh % Nodes % z(t)
    END DO
    CLOSE(10)
#endif
    
  END SUBROUTINE MarkSharpNodes
    

  
!------------------------------------------------------------------------------
!> Finds neighbours of the nodes in given direction.
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
    CALL Warn('FindNeigbourNodes','SIZE of Neighbours should equal Number of Nodes!')
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
  SUBROUTINE UpdateSolverMesh( Solver, Mesh, NoInterp )
!------------------------------------------------------------------------------
     TYPE( Mesh_t ), POINTER :: Mesh
     TYPE( Solver_t ), TARGET :: Solver
     LOGICAL, OPTIONAL :: NoInterp
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,n1,n2,DOFs
     LOGICAL :: Found, OptimizeBandwidth, GlobalBubbles, IsTransient
     TYPE(Matrix_t), POINTER   :: Matrix
     REAL(KIND=dp), POINTER :: Work(:)
     INTEGER, POINTER :: Permutation(:)
     TYPE(Variable_t), POINTER :: TimeVar, SaveVar, Var
     CHARACTER(:), ALLOCATABLE :: str
     LOGICAL :: DoInterp 
!------------------------------------------------------------------------------
     SaveVar => Solver % Variable
     DOFs = SaveVar % DOFs
!
!    Create matrix and variable structures for
!    current equation on the new mesh:
!    -----------------------------------------

     ! Backward compatibility
     DoInterp = .TRUE.
     IF(PRESENT(NoInterp)) THEN
       DoInterp = .NOT. NoInterp
     END IF

     Solver % Mesh => Mesh
     CALL SetCurrentMesh( CurrentModel, Mesh )

     IF  (DoInterp) THEN
       Solver % Variable => VariableGet( Mesh % Variables, &
           SaveVar % Name, ThisOnly = .FALSE. )
       CALL AllocateVector(Permutation, SIZE(Solver % Variable % Perm))
     ELSE
       ALLOCATE(Permutation(Mesh % NumberOfNodes + &
           Solver % Mesh % MaxEdgeDofs*Mesh % NumberOfEdges + &
           Solver % Mesh % MaxFaceDofs*Mesh % NumberOfFaces + &
           Solver % Mesh % MaxBDofs*Mesh % NumberOfBulkElements))
     END IF
     Permutation = 0
     
     GlobalBubbles = Solver % GlobalBubbles
     
     OptimizeBandwidth = ListGetLogical( Solver % Values, 'Optimize Bandwidth', Found )
     IF ( .NOT. Found ) OptimizeBandwidth = .TRUE.
     
     Matrix => CreateMatrix( CurrentModel, Solver, &
         Mesh, Permutation, DOFs, MATRIX_CRS, OptimizeBandwidth, &
         ListGetString( Solver % Values, 'Equation' ), &
         GlobalBubbles=GlobalBubbles)

     IF( ASSOCIATED( Matrix ) ) THEN
       Matrix % Symmetric = ListGetLogical( Solver % Values, &
           'Linear System Symmetric', Found )

       Matrix % Lumped = ListGetLogical( Solver % Values, &
           'Lumped Mass Matrix', Found )    
     END IF

     IF(.NOT. DoInterp) THEN
       Solver % Variable => VariableGet( Mesh % Variables, &
           SaveVar % Name, ThisOnly = .TRUE. )                     
       IF(.NOT. ASSOCIATED( Solver % Variable ) ) THEN
         CALL VariableAddVector( Mesh % Variables, Mesh, Solver, &
             SaveVar % Name, SaveVar % Dofs, Perm = Permutation )
         Solver % Variable => VariableGet( Mesh % Variables, &
             SaveVar % Name, ThisOnly = .TRUE. )                     
       END IF
         
       Solver % Variable % Perm => Permutation
       IF(.NOT. ASSOCIATED( Solver % Variable % perm) ) THEN
         CALL Fatal('UpdateSolverMesh','No Perm associated?!')
       END IF
       NULLIFY(Permutation)

       IsTransient = ( ListGetString( CurrentModel % Simulation,&
           'Simulation Type' ) == 'transient' ) 
       IF( IsTransient ) THEN
         n1 = SIZE( Solver % Variable % Values )
         IF ( Solver % TimeOrder == 2 ) THEN
           n2 = 7
         ELSE 
           n2 = MAX( Solver % Order, Solver % TimeOrder )
         END IF
         ALLOCATE( Solver % Variable % PrevValues(n1,n2) )
         Solver % Variable % PrevValues = 0.0_dp
       END IF         
     ELSE
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
       DEALLOCATE( Permutation )
     END IF

     Solver % Variable % Solver => Solver


     IF (ASSOCIATED(Matrix)) CALL AllocateVector( Matrix % RHS, Matrix % NumberOfRows )

     IF ( ASSOCIATED(SaveVar % EigenValues) ) THEN
       n = SIZE(SaveVar % EigenValues)

       IF ( n > 0 ) THEN
         Solver % NOFEigenValues = n
         CALL AllocateVector( Solver % Variable % EigenValues,n )
         CALL AllocateArray( Solver % Variable % EigenVectors, n, &
             SIZE(Solver % Variable % Values) ) 

         IF( Solver % Variable % Dofs > 1 ) THEN
           DO k=1,Solver % Variable % DOFs
             str = ComponentName( Solver % Variable % Name, k )
             Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )
             IF ( ASSOCIATED( Var ) ) THEN
               Var % EigenValues => Solver % Variable % EigenValues
               Var % EigenVectors =>  & 
                   Solver % Variable % EigenVectors(:,k::Solver % Variable % DOFs )
             END IF
           END DO
         END IF

         Solver % Variable % EigenValues  = 0.0d0
         Solver % Variable % EigenVectors = 0.0d0

         IF (ASSOCIATED(Matrix)) THEN
           CALL AllocateVector( Matrix % MassValues, SIZE(Matrix % Values) )
           Matrix % MassValues = 0.0d0
         END IF
       END IF
     ELSE IF ( ASSOCIATED( Solver % Matrix ) ) THEN
       IF( ASSOCIATED( Solver % Matrix % Force) ) THEN
         n1 = Matrix % NumberOFRows
         n2 = SIZE(Solver % Matrix % Force,2)
         ALLOCATE(Matrix % Force(n1,n2))
         Matrix % Force = 0.0d0
       END IF
     END IF

     IF (ASSOCIATED(Matrix)) THEN
       Solver % Matrix => Matrix
     ELSE
       NULLIFY(Solver % Matrix)
     END IF
     Solver % Mesh % Changed = .TRUE.

!------------------------------------------------------------------------------
  END SUBROUTINE UpdateSolverMesh
!------------------------------------------------------------------------------



  ! Create list of active elements for more speedy operation
  !-------------------------------------------------------------
  SUBROUTINE SetActiveElementsTable( Model, Solver, MaxDim, CreateInv )
    TYPE(Model_t)  :: Model
    TYPE(Solver_t) :: Solver
    INTEGER, OPTIONAL :: MaxDim
    LOGICAL, OPTIONAL :: CreateInv
    
    INTEGER :: i, n, Sweep, MeshDim 
    TYPE(Element_t), POINTER :: Element
    LOGICAL :: Found, HasFCT, Parallel
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(:), ALLOCATABLE :: EquationName
    
    IF( .NOT. ( Solver % Mesh % Changed .OR. Solver % NumberOfActiveElements <= 0 ) ) RETURN

    IF( ASSOCIATED( Solver % ActiveElements ) ) THEN
      DEALLOCATE( Solver % ActiveElements )
    END IF
    
    EquationName = ListGetString( Solver % Values, 'Equation', Found)
    IF( .NOT. Found ) THEN
      CALL Fatal('SetActiveElementsTable','Equation not present!')
    END IF

    CALL Info('SetActiveElementsTable',&
        'Creating active element table for: '//TRIM(EquationName),Level=12)

    HasFCT = ListGetLogical( Solver % Values, 'Linear System FCT', Found )

    Mesh => Solver % Mesh

    MeshDim = 0 
    Parallel = ( ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 

    
    DO Sweep = 0, 1    
      n = 0
      DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
        Element => Solver % Mesh % Elements(i)

        IF( Parallel ) THEN
          IF( .NOT.HasFCT .AND. Element % PartIndex /= ParEnv % myPE ) CYCLE
        END IF
          
        IF ( CheckElementEquation( Model, Element, EquationName ) ) THEN
          n = n + 1
          IF( Sweep == 0 ) THEN
            MeshDim = MAX( Element % TYPE % DIMENSION, MeshDim )
          ELSE
            Solver % ActiveElements(n) = i
          END IF
        END IF
      END DO
      
      IF( Sweep == 0 ) THEN
        Solver % NumberOfActiveElements = n
        IF( n == 0 ) EXIT
        ALLOCATE( Solver % ActiveElements( n ) )
      END IF
    END DO

    IF( n == 0 ) THEN
      CALL Info('SetActiveElementsTable','No active elements found',Level=12)    
      RETURN
    END IF
                
    IF( PRESENT( MaxDim ) ) MaxDim = MeshDim 

    IF( PRESENT( CreateInv ) ) THEN
      IF( CreateInv ) THEN
        CALL Info('SetActiveElementsTable','Creating inverse table for elemental variable permutation',Level=20)
        ALLOCATE( Solver % InvActiveElements( Mesh % NumberOfBulkElements &
            + Mesh % NumberOFBoundaryElements ) )

        Solver % InvActiveElements = 0
        DO i=1,Solver % NumberOfActiveElements
          Solver % InvActiveElements( Solver % ActiveElements(i) ) = i
        END DO
      END IF
    END IF
    
    CALL Info('SetActiveElementsTable','Number of active elements found : '//I2S(n),Level=12)    
    
  END SUBROUTINE SetActiveElementsTable



!------------------------------------------------------------------------------
!> Split a mesh equally to smaller pieces by performing a uniform split.
!> Also known as mesh multiplication. A 2D element splits into 4 elements of
!> same form, and 3D element into 8 elements. 
!> Currently works only for linear elements.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Split mesh with quadfaces into one with only triangle faces.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Sometimes we are lucky and the mesh includes similar elements that are
!> different only by their center point. If we then ensure that their local
!> numbering is the same we may use same finite element basis vectors for them.
!------------------------------------------------------------------------------
  SUBROUTINE SetEqualElementIndeces( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), TARGET :: Mesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: r0(:,:), r1(:,:)
    REAL(KIND=dp) :: eps, dist
    INTEGER, POINTER :: Indexes0(:)
    INTEGER, ALLOCATABLE :: Indexes1(:)
    INTEGER :: t,t0,i,j,n,n0,n1,na,nb,cnt(2)
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: SimilarElement(:)
    LOGICAL :: Similar
    CHARACTER(:), ALLOCATABLE :: str    
    
    n = Mesh % MaxElementNodes
    ALLOCATE(r0(n,3),r1(n,3),Indexes1(n))

    cnt = [1,0]
    na = Mesh % NumberOfBulkElements
    nb = Mesh % NumberOfBoundaryElements

    ALLOCATE(SimilarElement(na+nb))
    SimilarElement = 0
    
    DO t=1,na+nb
      Element => Mesh % Elements(t)
      Indexes0 => Element % NodeIndexes
      n = Element % Type % NumberOfNodes

      r1(1:n,1) = Mesh % Nodes % x(Indexes0)
      r1(1:n,2) = Mesh % Nodes % y(Indexes0)
      r1(1:n,3) = Mesh % Nodes % z(Indexes0)

      ! Compute distances from element center. 
      DO i=1,3
        r1(1:n,i) = r1(1:n,i) - SUM(r1(1:n,i))/n
      END DO

      ! Memorize the reference element. 
      IF(t==1) THEN
        r0 = r1
        n0 = n
        t0 = t
        eps = 1.0e-6 * SUM(ABS(r0))/n0
        CYCLE
      END IF

      Similar = .FALSE.
      IF(n == n0) THEN
        n1 = 0
        DO i=1,n
          DO j=1,n
            dist = SQRT(SUM((r1(j,:)-r0(i,:))**2))
            IF(dist < eps) THEN
              n1 = n1 + 1
              Indexes1(i) = Indexes0(j)
            END IF
          END DO
        END DO        
        IF(n1 == n) THEN
          Similar = .TRUE.
          cnt(1) = cnt(1) + 1
          IF(ANY(Indexes0(1:n) /= Indexes1(1:n))) THEN
            cnt(2) = cnt(2) + 1
            Indexes0(1:n) = Indexes1(1:n)
          END IF
          SimilarElement(t) = t0
        END IF
      END IF

      ! Create new reference!
      IF(.NOT. Similar) THEN
        r0 = r1
        n0 = n
        t0 = t
      END IF

      IF( t == na .OR. t == na + nb ) THEN
        IF( t == na ) THEN
          str = 'bulk'
          n = na
        ELSE
          str = 'boundary'
          n = nb
        END IF
        CALL Info('SetEqualElementIndeces','Number of Similar '//TRIM(str)//' elements '&
            //I2S(cnt(1))//' (out of '//I2S(n)//')')
        CALL Info('SetEqualElementIndeces','Number of altered '//TRIM(str)//' elements')
        cnt = 0
      END IF
    END DO

    DEALLOCATE( SimilarElement ) 
    
  END SUBROUTINE SetEqualElementIndeces
    
  


!------------------------------------------------------------------------------
  SUBROUTINE SetCurrentMesh( Model, Mesh )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t),  POINTER :: Mesh
!------------------------------------------------------------------------------

    IF(.NOT. ASSOCIATED(Mesh) ) THEN
      CALL Fatal('SetCurrentMesh','Target mesh is not associated!')
    END IF

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


!----------------------------------------------------------------------------------
  SUBROUTINE DisplaceMesh( Mesh, Update, sgn, Perm, DOFs, StabRecomp, UpdateDirs )
!----------------------------------------------------------------------------------
    TYPE(Mesh_t) , POINTER :: Mesh 
    REAL(KIND=dp) :: Update(:)
    INTEGER :: DOFs,sgn,Perm(:)
    LOGICAL, OPTIONAL :: StabRecomp
    INTEGER, OPTIONAL :: UpdateDirs

    INTEGER :: i,k,dim
    LOGICAL :: StabFlag

    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element

    IF ( PRESENT( UpdateDirs ) ) THEN
      dim = UpdateDirs
    ELSE
      dim = DOFs
    END IF

    DO i=1,MIN( SIZE(Perm), SIZE(Mesh % Nodes % x) )
       k = Perm(i)
       IF ( k > 0 ) THEN
         k = DOFs * (k-1)
         Mesh % Nodes % x(i)   = Mesh % Nodes % x(i) + sgn * Update(k+1)
         IF ( dim > 1 ) &
           Mesh % Nodes % y(i) = Mesh % Nodes % y(i) + sgn * Update(k+2)
         IF ( dim > 2 ) &
           Mesh % Nodes % z(i) = Mesh % Nodes % z(i) + sgn * Update(k+3)
        END IF
    END DO

    StabFlag = .TRUE.
    IF ( PRESENT( StabRecomp ) ) StabFlag = StabRecomp

    IF ( sgn == 1 .AND. StabFlag ) THEN
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
       ! Check if last node matches global max node
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
!> Assign local number of edge to given boundary element. Also copies all 
!> p element attributes from element edge to boundary edge.
!------------------------------------------------------------------------------
  SUBROUTINE AssignLocalNumber( EdgeElement, Element, Mesh, NoPE )
!------------------------------------------------------------------------------
    USE PElementMaps, ONLY : getFaceEdgeMap 
    IMPLICIT NONE

    ! Parameters
    TYPE(Mesh_t) :: Mesh                     !< Finite element mesh containing faces and edges.
    TYPE(Element_t), POINTER :: EdgeElement  !< Edge element to which assign local number
    TYPE(Element_t), POINTER :: Element      !< Bulk element with some global numbering to use to assign local number
    LOGICAL, OPTIONAL :: NoPE
!------------------------------------------------------------------------------
    ! Local variables

    INTEGER i,j,k,n,edgeNumber, numEdges, bMap(4), bIndex(4)
    TYPE(Element_t), POINTER :: Edge
    LOGICAL :: EvalPE

    EvalPE = .TRUE.
    IF(PRESENT(NoPE)) EvalPE = .NOT.NoPE
    
    ! Get number of points, edges or faces
    numEdges = 0
    SELECT CASE (Element % TYPE % DIMENSION)
    CASE (0,1)
      RETURN
    CASE (2)
       numEdges = Element % TYPE % NumberOfEdges
    CASE (3)   
       numEdges = Element % TYPE % NumberOfFaces
    CASE DEFAULT
      CALL Fatal('AssignLocalNumber','Unsupported Element dim: '//I2S(Element % TYPE % DIMENSION))
    END SELECT

    ! If edges have not been created, stop search. This should not happen, actually.
    IF (.NOT. ASSOCIATED(Element % EdgeIndexes)) THEN
      CALL Warn('AssignLocalNumber','Edge indexes for element not associated!')
      RETURN
    END IF
        
    ! For each edge or face in element try to find local number
    DO edgeNumber=1, numEdges
      Edge => GetElementEntity(Element,edgeNumber,Mesh)
      
      ! Edge element not found. This should not be possible, unless there
      ! is an error in the mesh read in process..
      IF (.NOT. ASSOCIATED(Edge)) THEN
        CALL Fatal('MeshUtils::AssignLocalNumber','Edge element not found')
      END IF
      
      n = 0
      ! For each element node
      DO i=1, Edge % TYPE % NumberOfNodes
        ! For each node in edge element
        DO j=1, EdgeElement % TYPE % NumberOfNodes
          ! If edge and edgeelement node match increment counter
          IF (Edge % NodeIndexes(i) == EdgeElement % NodeIndexes(j)) THEN
            n = n + 1
            EXIT
          END IF
        END DO
      END DO

      ! If all nodes are on boundary, edge/face was found
      IF (n == EdgeElement % TYPE % NumberOfNodes) THEN
        IF(EvalPE) THEN
          EdgeElement % PDefs % localNumber = edgeNumber
          EdgeElement % PDefs % LocalParent => Element
        END IF

        ! Change ordering of global nodes to match that of element
        bMap = getElementBoundaryMap( Element, edgeNumber )
        Bindex(1:n) = Element % NodeIndexes(bMap(1:n))

        k = 0
        DO i=1, Edge % TYPE % NumberOfNodes
          DO j=1, EdgeElement % TYPE % NumberOfNodes
            IF (Edge % NodeIndexes(i) == bIndex(j)) THEN
              k = k + 1
              EXIT
            END IF
          END DO
        END DO
        
        ! Ok, reorder the nodal to comply with the mapping.
        ! Do not do this if we would not just reorder but also loose some nodes!
        IF(k==n) THEN
          EdgeElement % NodeIndexes(1:n) = Bindex(1:n)
        ELSE
#if 0
          PRINT *,'Element Types: ',Element % TYPE % ElementCode, EdgeElement % TYPE % ElementCode, numEdges
          IF(ASSOCIATED(Element % Pdefs)) PRINT *,'Element TetraType:',Element % PDefs % TetraType 
          PRINT *,'Element:',Element % NodeIndexes
          PRINT *,'EdgeA:  ',EdgeElement % NodeIndexes
          PRINT *,'EdgeB:  ',Edge % NodeIndexes
          PRINT *,'Element hits:',EvalPE, n,k
          PRINT *,'BoundaryMap:',bmap(1:n)
#endif
          CALL Warn('AssignLocalNumber','Skipped mapping as we would loose nodes!')
        END IF

        ! Copy misc attributes of edge element to boundary element
        IF(EvalPE) THEN
          EdgeElement % PDefs % isEdge = Edge % PDefs % isEdge
          EdgeElement % PDefs % GaussPoints = Edge % PDefs % GaussPoints
          EdgeElement % PDefs % P = Edge % PDefs % P
        END IF

        !(and boundary bubble dofs)
        EdgeElement % BDOFs = MAX(EdgeElement % BDOFs, Edge % BDOFs)

        
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
    IF(.NOT.ASSOCIATED(EdgeElement % PDefs % LocalParent)) THEN
      CALL Warn('MeshUtils::AssignLocalNumber','Unable to find local edge '//I2S(EdgeElement % ElementIndex))
    END IF

        
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
        CALL Fatal('AssignLocalNumber::GetElementEntity',&
            'Impossible Element dim: '//I2S(Element % TYPE % DIMENSION))
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
      OptimizeBW, Perm, LocalNodes, MaskOnBulk, RequireLogical, &
      ParallelComm, BreakLoop )
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
    LOGICAL, OPTIONAL :: ParallelComm
    LOGICAL, OPTIONAL :: BreakLoop
!------------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm(:), Neighbours(:)
    INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:)
    TYPE(ListMatrix_t), POINTER :: ListMatrix(:)
    INTEGER :: t,i,j,k,l,m,k1,k2,n,p,q,e1,e2,f1,f2,This,bf_id,nn,t0,ii(ParEnv % PEs)
    INTEGER :: ierr, status(MPI_STATUS_SIZE), NewDofs
    LOGICAL :: Flag, Found, FirstRound, MaskIsLogical, Hit, Parallel
    LOGICAL, ALLOCATABLE :: IsNeighbour(:)
    INTEGER :: Indexes(30), ElemStart, ElemFin, Width, BreakNode
    TYPE(ListMatrixEntry_t), POINTER :: CList, Lptr
    TYPE(Element_t), POINTER :: CurrentElement,Elm
    REAL(KIND=dp) :: MinDist, Dist
!------------------------------------------------------------------------------

    IF(PRESENT(ParallelComm)) THEN
      Parallel = ParallelComm .AND. ( ParEnv % PEs > 1 )
    ELSE
      Parallel = ParEnv % PEs > 1
    END IF

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

    IF( ElemFin - ElemStart <= 0 .AND. .NOT. Parallel) THEN
       LocalNodes = 0
       RETURN
    END IF

    k = 0
    Perm = 0
    FirstRound = .TRUE.
    BreakNode = 0
    t0 = 0
    
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
       
       n = CurrentElement % TYPE % NumberOfNodes
       Indexes(1:n) = CurrentElement % NodeIndexes(1:n)
       
       IF( FirstRound ) THEN
         ! Just plainly create the permutation
         DO i=1,n
             j = Indexes(i)
             IF ( Perm(j) == 0 ) THEN
                k = k + 1
                Perm(j) = k
             END IF
          END DO
        ELSE
          ! Create the list matrix for the connectivity in order to minimize the bandwidth
          DO i=1,n
             k1 = Perm(Indexes(i))
             IF ( k1 <= 0 ) CYCLE
             DO j=1,n
                k2 = Perm(Indexes(j))
                IF ( k2 <= 0 ) CYCLE
                IF( k1 == BreakNode .OR. k2 == BreakNode ) THEN
                  IF( t0 == 0 ) t0 = t
                  IF( t0 /= t ) THEN
                    PRINT *,'breaking connection between:',k1,k2
                    CYCLE
                  END IF
                END IF
                Lptr => List_GetMatrixIndex( ListMatrix,k1,k2 )
             END DO
          END DO
       END IF
    END DO
    LocalNodes = k

    ! In parallel case, detect nodes which are shared with another partition
    ! which may not have an element on this boundary
    ! Code borrowed from CommunicateLinearSystemTag
    !------------------------------------------------------------------------------
    IF( Parallel ) THEN

      ALLOCATE( IsNeighbour(ParEnv % PEs), fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )

      nn = MeshNeighbours(Mesh, IsNeighbour)
      nn = 0
      ineigh = 0
      DO i=0, ParEnv % PEs-1
        k = i+1
        IF(i==ParEnv % myPE) CYCLE
        IF(.NOT. IsNeighbour(k) ) CYCLE
        nn = nn + 1
        fneigh(nn) = k
        ineigh(k) = nn
      END DO

      n = COUNT(Perm > 0 .AND. Mesh % ParallelInfo % GInterface)
      ALLOCATE( s_e(n, nn ), r_e(n) )

      CALL CheckBuffer( nn*3*n )

      ii = 0
      DO i=1, Mesh % NumberOfNodes
        IF(Perm(i) > 0 .AND. Mesh % ParallelInfo % GInterface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k == ParEnv % MyPE ) CYCLE
            k = k + 1
            k = ineigh(k)
            IF ( k> 0) THEN
              ii(k) = ii(k) + 1
              s_e(ii(k),k) = Mesh % ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
        END IF
      END DO

      DO i=1, nn
        j = fneigh(i)
        CALL MPI_BSEND( ii(i),1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD,ierr )
        IF( ii(i) > 0 ) THEN
          CALL MPI_BSEND( s_e(1:ii(i),i),ii(i),MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,ierr )
        END IF
      END DO

      NewDofs = 0

      DO i=1, nn
        j = fneigh(i)
        CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
        IF ( n>0 ) THEN
          IF( n>SIZE(r_e)) THEN
            DEALLOCATE(r_e)
            ALLOCATE(r_e(n))
          END IF

          CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
          DO j=1,n
            k = SearchNode( Mesh % ParallelInfo, r_e(j), Order=Mesh % ParallelInfo % Gorder )
            IF ( k>0 ) THEN
              IF(.NOT. Perm(k) > 0) THEN
                NewDofs = NewDofs + 1
                Perm(k) = LocalNodes + NewDofs
              END IF
            END IF
          END DO
        END IF
      END DO
      DEALLOCATE(s_e, r_e )

      LocalNodes = LocalNodes + NewDofs
    END IF

    ! Don't optimize bandwidth for parallel cases
    IF( Parallel .OR. .NOT. OptimizeBW ) RETURN

    IF(FirstRound) THEN
       ! Allocate space 
       NULLIFY( ListMatrix )
       ListMatrix => List_AllocateMatrix(LocalNodes)
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

       ! Minimizing the bandwidth of a closed loop is impossible.
       ! So let us break the loop on one node. 
       IF(PRESENT(BreakLoop)) THEN
         IF(BreakLoop) BreakNode = 1
       END IF
       
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
    REAL(kind=dp) :: BoundingBox(6), eps2, eps1 = 1d-3, GlobalEps, LocalEps
    CHARACTER(:), ALLOCATABLE :: MaskName


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
      IF(.NOT. Stat ) IsRecursive = .TRUE.

      LocalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Local Epsilon', Stat )
      IF(.NOT. stat) LocalEps = 1.0d-10

      GlobalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Global Epsilon', Stat ) 
      IF(.NOT. stat) THEN
        IF( IsRecursive ) THEN
          GlobalEps = 2.0d-10
        ELSE
          GlobalEps = 1.0d-4
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


  !> Calculate the number of separature pieces in a serial mesh.
  !> This could be used to detect problems in mesh when suspecting
  !> floating parts not fixed by any BC, for example.
  !---------------------------------------------------------------------------------
  SUBROUTINE CalculateMeshPieces( Mesh, ElementMode, PieceIndex)

    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: ElementMode
    INTEGER, OPTIONAL :: PieceIndex(:)

    LOGICAL :: Ready
    INTEGER :: i,j,k,n,t,t2,k2,MinIndex,MaxIndex,Loop,NoPieces
    INTEGER, ALLOCATABLE :: MeshPiece(:),PiecePerm(:)
    TYPE(Element_t), POINTER :: Element, Element2
    INTEGER, POINTER :: Indexes(:)
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Faces(:)
    LOGICAL :: ElemMode, Found
    
    IF( ParEnv % PEs > 1 ) THEN
      CALL Warn('CalculateMeshPieces','Implemented only for serial meshes!')
    END IF

    ElemMode = .FALSE.
    IF( PRESENT(ElementMode) ) THEN
      ElemMode = ElementMode
    END IF

    IF( ElemMode ) THEN
      n = Mesh % NumberOfBulkElements
    ELSE   
      n = Mesh % NumberOfNodes
    END IF
    ALLOCATE( MeshPiece( n ) ) 
    MeshPiece = 0

    ! Only set the piece for the nodes that are used by some element
    ! For others the marker will remain zero. 
    DO t = 1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)        
      IF( ElemMode ) THEN
        MeshPiece( t ) = 1
      ELSE      
        Indexes => Element % NodeIndexes
        MeshPiece( Indexes ) = 1
      END IF
    END DO
    j = 0
    DO i = 1, n
      IF( MeshPiece(i) > 0 ) THEN
        j = j + 1
        MeshPiece(i) = j
      END IF
    END DO

    IF(n>j) THEN
      CALL Info('CalculateMeshPieces',&
          'Number of non-body nodes in mesh is '//I2S(n-j),Level=5)
    END IF
      
    ! We go through the elements and set all the piece indexes to minimimum index
    ! until the mesh is unchanged. Thereafter the whole piece will have the minimum index
    ! of the piece.
    Ready = .FALSE.
    Loop = 0
    DO WHILE(.NOT. Ready) 
      Ready = .TRUE.
      DO t = 1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)        
        
        IF( ElemMode ) THEN
          k = MeshPiece(t)
          IF( Mesh % MeshDim == 2 ) THEN
            DO i=1, Element % TYPE % NumberOfEdges
              DO j=1,2
                IF(j==1) THEN
                  Element2 => Mesh % Edges(Element % EdgeIndexes(i)) % BoundaryInfo % Left
                ELSE
                  Element2 => Mesh % Edges(Element % EdgeIndexes(i)) % BoundaryInfo % Right
                END IF
                IF(.NOT. ASSOCIATED(Element2) ) CYCLE
                t2 = Element2 % ElementIndex
                IF(t==t2) CYCLE
                k2 = MeshPiece(t2)
                IF(k2 /= k ) THEN
                  Ready = .FALSE.
                  IF( k2 < k ) THEN
                    k = k2 
                    MeshPiece(t) = k2
                  ELSE
                    MeshPiece(t2) = k
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO i=1, Element % TYPE % NumberOfFaces
              DO j=1,2
                IF(j==1) THEN
                  Element2 => Mesh % Faces(Element % FaceIndexes(i)) % BoundaryInfo % Left
                ELSE
                  Element2 => Mesh % Faces(Element % FaceIndexes(i)) % BoundaryInfo % Right
                END IF
                IF(.NOT. ASSOCIATED(Element2) ) CYCLE
                t2 = Element2 % ElementIndex
                IF(t==t2) CYCLE
                k2 = MeshPiece(t2)
                IF(k2 /= k ) THEN
                  Ready = .FALSE.
                  IF( k2 < k ) THEN
                    k = k2 
                    MeshPiece(t) = k2
                  ELSE
                    MeshPiece(t2) = k
                  END IF
                END IF
              END DO
            END DO
          END IF
        ELSE
          Indexes => Element % NodeIndexes          
          MinIndex = MINVAL( MeshPiece( Indexes ) )
          MaxIndex = MAXVAL( MeshPiece( Indexes ) )
          IF( MaxIndex > MinIndex ) THEN
            MeshPiece( Indexes ) = MinIndex
            Ready = .FALSE.
          END IF
        END IF
      END DO
      Loop = Loop + 1
    END DO
    CALL Info('CalculateMeshPieces','Mesh coloring loops: '//I2S(Loop),Level=6)

    ! Compute the true number of different pieces
    MaxIndex = MAXVAL( MeshPiece )
    IF( MaxIndex == 1 ) THEN
      NoPieces = 1
      IF(PRESENT(PieceIndex)) PieceIndex = 1
    ELSE
      ALLOCATE( PiecePerm( MaxIndex ) ) 
      PiecePerm = 0
      NoPieces = 0
      DO i = 1, n
        j = MeshPiece(i) 
        IF( j == 0 ) CYCLE
        IF( PiecePerm(j) == 0 ) THEN
          NoPieces = NoPieces + 1
          PiecePerm(j) = NoPieces 
        END IF
      END DO
      ! Use the compact numbering of mesh pieces
      DO i=1,n
        j = MeshPiece(i)
        IF(j>0) MeshPiece(i) = PiecePerm(j)
      END DO
      IF(PRESENT(PieceIndex)) PieceIndex = MeshPiece
    END IF
    CALL Info('CalculateMeshPieces',&
        'Number of separate pieces in mesh is '//I2S(NoPieces),Level=5)
    
    IF(PRESENT(PieceIndex)) RETURN
    
    i = ListGetInteger( CurrentModel % Simulation,'Desired Mesh Pieces',Found )
    IF( Found ) THEN
      IF( i == NoPieces ) THEN
        CALL Info('CalculateMeshPieces','Number of pieces agree with the requested '//I2S(i))
        RETURN
      ELSE
        CALL Fatal('CalculateMeshPieces','Number of pieces differ from the requested '//I2S(i))
      END IF
    END IF

    ! No point to create piece of just ones
    IF( NoPieces == 1 ) RETURN
    
    ! Save the mesh piece field to > mesh piece < 
    Var => VariableGet( Mesh % Variables,'Mesh Piece' )
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      IF( ElemMode ) THEN
        CALL VariableAddVector ( Mesh % Variables,Mesh, CurrentModel % Solver,'Mesh Piece', &
            VarType = Variable_on_elements )
      ELSE
        CALL VariableAddVector ( Mesh % Variables,Mesh, CurrentModel % Solver,'Mesh Piece' )
      END IF
      Var => VariableGet( Mesh % Variables,'Mesh Piece' )
    END IF

    IF( .NOT. ASSOCIATED( Var ) ) THEN
      CALL Fatal('CalculateMeshPieces','Could not get handle to variable > Mesh Piece <')
    END IF

    DO i = 1, n
      j = i
      IF( ASSOCIATED( Var % Perm ) ) THEN
        j = Var % Perm( i ) 
        IF( j == 0 ) CYCLE
      END IF
      Var % Values( j ) = 1.0_dp * MeshPiece( i ) 
    END DO
    CALL Info('CalculateMeshPieces','Creating variable showing the non-connected domains: mesh piece',Level=5)
  
  END SUBROUTINE CalculateMeshPieces
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Compute radius of rotor using only topology information.
!> Assumes that axis of rotation is z-axis. 
!------------------------------------------------------------------------------
  FUNCTION DetermineRotorRadius(Mesh) RESULT( Radius )
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Radius
    
    INTEGER, ALLOCATABLE :: PieceIndex(:)
    INTEGER :: i,imin,n
    REAL(KIND=dp) :: r2,rmin,rmax

    Radius = -1.0_dp
    n = Mesh % NumberOfNodes
    ALLOCATE(PieceIndex(n))
    PieceIndex = 0
    CALL CalculateMeshPieces( Mesh, PieceIndex = PieceIndex )
    IF( MAXVAL(PieceIndex) /= 2) RETURN

    ! Find minimum radius nodes i.e. center node
    rmin = HUGE(rmin)
    imin = 0
    DO i=1,n
      r2 = Mesh % Nodes % x(i)**2 + Mesh % Nodes % y(i)**2
      IF(r2<rmin) THEN
        rmin = r2
        imin = i
      END IF
    END DO

    ! Find the maximum radius in the same piece i.e. rotor radius
    rmax = 0.0_dp
    DO i=1,n
      IF(PieceIndex(i) /= PieceIndex(imin)) CYCLE
      r2 = Mesh % Nodes % x(i)**2 + Mesh % Nodes % y(i)**2
      rmax = MAX(rmax,r2)
    END DO
    Radius = SQRT(rmax)             
    
  END FUNCTION DetermineRotorRadius
!------------------------------------------------------------------------------
  

!--------------------------------------------------------------------------
!> This subroutine finds the structure of an extruded mesh even though it is 
!> given in an unstructured format. The routine may be used by some special
!> solvers that employ the special character of the mesh.
!> The extrusion is found for a given direction and for each node the corresponding 
!> up and down, and thereafter top and bottom node is computed.
!-----------------------------------------------------------------------------
!---------------------------------------------------------------


 
  !> Find the node closest to the given coordinate. 
  !> The linear search only makes sense for a small number of points. 
  !> Users include saving routines of pointwise information. 
  !-----------------------------------------------------------------
  FUNCTION ClosestNodeInMesh(Mesh,Coord,MinDist,DoParallel) RESULT ( NodeIndx )
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coord(3)
    REAL(KIND=dp), OPTIONAL :: MinDist
    LOGICAL, OPTIONAL :: DoParallel
    INTEGER :: NodeIndx

    REAL(KIND=dp) :: Dist2,MinDist2,ParDist2, NodeCoord(3)
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
    
    ! In parallel only return a hit in the correct partition.
    IF(PRESENT(DoParallel)) THEN
      IF( DoParallel ) THEN
        ParDist2 = ParallelReduction(MinDist2,1)
        IF(ABS(ParDist2-MinDist2) > 1.0e-20 ) THEN
          NodeIndx = 0
        END IF
      END IF
    END IF
      
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
        ! since afterwards there is no way of deciding which one was closer.
        !--------------------------------------------------------------------------------------
        IF( ParallelCands > 1.5_dp ) THEN
          Hit = PointInElement( Element, ElementNodes, &
              Coords, LocalCoords, GlobalEps = 1.0d-3, LocalEps=1.0d-4 )
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
    CHARACTER(:), ALLOCATABLE :: Method
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


      ! Choose the nodes within the cones in the given three directions
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
        ig => Mesh % ParallelInfo % GInterface
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

        ! Do not measure distance to unset nodes!
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
    INTEGER :: i, j, k, n, NoNodes, NoElements, ActiveDirection, Order, BodyId, ne
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t),POINTER :: elmt
    REAL(KIND=dp) :: MeshVector(3), Length, Coord(3)
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    CHARACTER(:), ALLOCATABLE :: MeshName
    
!------------------------------------------------------------------------------
    Mesh => NULL()
    IF ( .NOT. ASSOCIATED( Params ) ) RETURN
    Mesh => AllocateMesh()

    CALL Info('CreateLineMesh','Creating 1D mesh on-the-fly')

!   Read in the parameters defining a uniform 1D mesh
!--------------------------------------------------------------    
    Order = ListGetInteger( Params,'1D Element Order',Found,minv=1,maxv=2)
    NoElements = ListGetInteger( Params,'1D Number Of Elements',minv=1)
    Length = ListGetConstReal( Params,'1D Mesh Length',Found)
    IF(.NOT. Found) Length = 1.0_dp
    ActiveDirection = ListGetInteger( Params,'1D Active Direction',Found,minv=-3,maxv=3)
    IF(.NOT.Found) ActiveDirection = 1
    BodyId = ListGetInteger( Params,'1D Body Id',Found,minv=1)
    IF(.NOT. Found) BodyId = 1
    MeshName = ListGetString( Params,'1D Mesh Name',Found)
    IF(.NOT. Found) MeshName = '1d_mesh'
    
    Mesh % Name = TRIM(MeshName)
    Mesh % OutputActive = .FALSE.

!   Compute the resulting mesh parameters
!--------------------------------------------------------------
    ne = Order + 1
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

    ALLOCATE( w(0:NoNodes-1) )
    
    CALL UnitSegmentDivision( w, NoNodes-1, Params )
    
    DO i=1, NoNodes
      Coord = MeshVector * w(i-1)

      x(i) = Coord(1)
      y(i) = Coord(2)
      z(i) = Coord(3)
    END DO
    

!   Define elements
!   -------------------------------
    CALL AllocateVector( Mesh % Elements, NoElements )

    Elmt => GetElementType( 200 + ne )

    DO i=1,NoElements
      Element => Mesh % Elements(i)      
      Element % TYPE => Elmt
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()     
      Element % ElementIndex = i

      CALL AllocateVector( Element % NodeIndexes, ne )
      Element % Ndofs = ne ! TO DO: This is not consistent for "Element = n:N", with N>1

      Element % NodeIndexes(1) = (i-1)*Order + 1
      Element % NodeIndexes(2) = i*Order + 1

      DO j=3,ne
        Element % NodeIndexes(j) = (i-1)*Order + j-1
      END DO
      
      Element % BodyId = BodyId
      Element % PartIndex = ParEnv % myPE
    END DO
    
!   Update new mesh node count:
!   ---------------------------

    Mesh % NumberOfNodes = NoNodes
    Mesh % Nodes % NumberOfNodes = NoNodes
    Mesh % NumberOfBulkElements = NoElements
    Mesh % MaxElementNodes = ne
    Mesh % MaxElementDOFs = ne
    Mesh % MeshDim = 1

    CALL SetMeshMaxDOFs(Mesh)

    
    WRITE(Message,'(A,I0)') 'Number of elements created: ',NoElements
    CALL Info('CreateLineMesh',Message)

    WRITE(Message,'(A,I0)') 'Number of nodes created: ',NoNodes
    CALL Info('CreateLineMesh',Message)
 
    CALL Info('CreateLineMesh','All done',Level=20)

  END FUNCTION CreateLineMesh

  !Creates a regular 2D mesh of 404 elements
  !The resulting mesh has no boundary elements etc for now
  !Should only be used for e.g. mesh to mesh interpolation
  FUNCTION CreateRectangularMesh(Params) RESULT(Mesh)

!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    REAL(KIND=dp) :: min_x, max_x, min_y, max_y, dx, dy
    INTEGER :: i, j, k, n, counter, nnx, nny, nex, ney, &
         NoNodes, NoElements, col, row
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t),POINTER :: elmt
    REAL(KIND=dp) :: MeshVector(3), Length, Coord(3)
    CHARACTER(*), PARAMETER :: FuncName="CreateRectangularMesh"

!------------------------------------------------------------------------------
    Mesh => NULL()
    IF ( .NOT. ASSOCIATED( Params ) ) RETURN
    Mesh => AllocateMesh()

    CALL Info(FuncName,'Creating 2D mesh on-the-fly')

    !Get parameters from valuelist
    min_x = ListGetConstReal(Params, "Grid Mesh Min X",UnfoundFatal=.TRUE.)
    max_x = ListGetConstReal(Params, "Grid Mesh Max X",UnfoundFatal=.TRUE.)
    min_y = ListGetConstReal(Params, "Grid Mesh Min Y",UnfoundFatal=.TRUE.)
    max_y = ListGetConstReal(Params, "Grid Mesh Max Y",UnfoundFatal=.TRUE.)
    dx    = ListGetConstReal(Params, "Grid Mesh dx",UnfoundFatal=.TRUE.)
    dy    = ListGetConstReal(Params, "Grid Mesh dy",Found)
    IF(.NOT. Found) dy = dx

    IF(max_x <= min_x .OR. max_y <= min_y .OR. dx <= 0.0_dp .OR. dy <= 0.0_dp) &
         CALL Fatal(FuncName, "Bad Grid Mesh parameters!")

    !number of nodes in x and y direction (and total)
    nnx = FLOOR((max_x - min_x) / dx) + 1
    nny = FLOOR((max_y - min_y) / dy) + 1
    NoNodes = nnx * nny

    !number of elements in x and y direction (and total)
    nex = nnx - 1
    ney = nny - 1
    NoElements = nex * ney


!   Define nodal coordinates
!   -------------------------------
    CALL AllocateVector( Mesh % Nodes % x, NoNodes )
    CALL AllocateVector( Mesh % Nodes % y, NoNodes )
    CALL AllocateVector( Mesh % Nodes % z, NoNodes )
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z

    z = 0.0_dp !2D

    !Define node positions
    counter = 0
    DO i=1,nnx
      DO j=1,nny
        counter = counter + 1
        x(counter) = min_x + (i-1)*dx
        y(counter) = min_y + (j-1)*dy
      END DO
    END DO

!   Define elements
!   -------------------------------
    CALL AllocateVector( Mesh % Elements, NoElements )

    Elmt => GetElementType( 404 )

    DO i=1,NoElements
      Element => Mesh % Elements(i)
      Element % TYPE => Elmt
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()
      Element % ElementIndex = i
      CALL AllocateVector( Element % NodeIndexes, 4 )
      Element % Ndofs = 4 ! TO DO: This is not consistent for "Element = n:N", with N>1

      col = MOD(i-1,ney)
      row = (i-1)/ney

      !THIS HERE NEEDS FIXED!!!!!
      Element % NodeIndexes(1) = (row * nny) + col + 1
      Element % NodeIndexes(2) = (row * nny) + col + 2
      Element % NodeIndexes(4) = ((row+1) * nny) + col + 1
      Element % NodeIndexes(3) = ((row+1) * nny) + col + 2

      Element % BodyId = 1
      Element % PartIndex = ParEnv % myPE
    END DO

!   Update new mesh node count:
!   ---------------------------

    Mesh % NumberOfNodes = NoNodes
    Mesh % Nodes % NumberOfNodes = NoNodes
    Mesh % NumberOfBulkElements = NoElements
    Mesh % MaxElementNodes = 4
    Mesh % MaxElementDOFs = 4
    Mesh % MeshDim = 2

  END FUNCTION CreateRectangularMesh


  !> Calculate body average for a discontinuous galerkin field.
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
    INTEGER :: n,i,j,k,l,nodeind,dgind, Nneighbours
    REAL(KIND=dp) :: AveHits
    LOGICAL, ALLOCATABLE :: IsNeighbour(:)
    LOGICAL :: Parallel

    
    IF(.NOT. ASSOCIATED(var)) RETURN
    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) RETURN

    IF( Var % DgAveraged ) THEN
      CALL Info('CalculateBodyAverage','Nodal average already computed for: '&
          //TRIM(Var % Name), Level=15)
      RETURN
    END IF
    
    IF( BodySum ) THEN
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal sum for: '&
          //TRIM(Var % Name), Level=8)
    ELSE
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal average for: '&
          //TRIM(Var % Name), Level=8)
    END IF

    Parallel = (ParEnv % PEs > 1 ) .AND. ( .NOT. Mesh % SingleMesh ) 
    
    
    n = Mesh % NumberOfNodes
    ALLOCATE( BodyCount(n), BodyAverage(n), IsNeighbour(Parenv % PEs) )
  
    
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
          ! This is just low priority info on the averaging
          IF( InfoActive(25) ) THEN
            j = COUNT(BodyCount > 0) 
            IF( j > 0 ) THEN
              AveHits = 1.0_dp * SUM( BodyCount ) / j
              WRITE(Message,'(A,ES12.3)') 'In body '//I2S(i)//' average hit count is: ',AveHits
              CALL Info('CalculateBodyAverage',Message) 
              WRITE(Message,'(A,2I0)') 'In body '//I2S(i)//' hit count range is: ',&
                  MINVAL(BodyCount,BodyCount>0), MAXVAL(BodyCount)
              CALL Info('CalculateBodyAverage',Message) 
            END IF
          END IF
        END IF
          
        IF( Parallel ) THEN
          Nneighbours = MeshNeighbours(Mesh, IsNeighbour)
          CALL SendInterface(); CALL RecvInterface()
        END IF

        j = COUNT( BodyCount > 0 )
        IF( j == 0 ) CYCLE
        
        ! Do not average weighted quantities (like nodal forces) - they should only be summed.
        ! But do average all other quantities. 
        IF( .NOT. BodySum ) THEN
          DO j=1,n
            IF( BodyCount(j) > 0 ) BodyAverage(j) = BodyAverage(j) / BodyCount(j)
          END DO
        END IF

        ! Now copy the average values to the DG field
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

    Var % DgAveraged = .TRUE.
    
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
         IF(.NOT.Mesh % ParallelInfo % GInterface(i)) CYCLE
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
         IF(.NOT.Mesh % ParallelInfo % GInterface(i)) CYCLE
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
         IF(.NOT. isNeighbour(i)) CYCLE

         CALL MPI_BSEND( cnt(i),1,MPI_INTEGER,i-1,1310,ELMER_COMM_WORLD,ierr )
         IF(cnt(i)>0) THEN
           CALL MPI_BSEND( buf(i) % gdof,cnt(i),MPI_INTEGER,i-1,1311,ELMER_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % ival,cnt(i),MPI_INTEGER,i-1,1312,ELMER_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % dval,cnt(i),MPI_DOUBLE_PRECISION,i-1,1313,ELMER_COMM_WORLD,ierr )
         END IF
       END DO
     END SUBROUTINE SendInterface


     SUBROUTINE RecvInterface()
       INTEGER, ALLOCATABLE :: gdof(:), ival(:)
       REAL(KIND=dp), ALLOCATABLE :: dval(:)
       INTEGER :: i,j,k,ierr, cnt, status(MPI_STATUS_SIZE)

       DO i=1,ParEnv % PEs

         IF(.NOT.isNeighbour(i)) CYCLE

         CALL MPI_RECV( cnt,1,MPI_INTEGER,i-1,1310,ELMER_COMM_WORLD,status,ierr )
         IF(cnt>0) THEN
           ALLOCATE( gdof(cnt), ival(cnt), dval(cnt) )
           CALL MPI_RECV( gdof,cnt,MPI_INTEGER,i-1,1311,ELMER_COMM_WORLD,status,ierr )
           CALL MPI_RECV( ival,cnt,MPI_INTEGER,i-1,1312,ELMER_COMM_WORLD,status,ierr )
           CALL MPI_RECV( dval,cnt,MPI_DOUBLE_PRECISION,i-1,1313,ELMER_COMM_WORLD,status,ierr )

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
       CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
     END SUBROUTINE RecvInterface

  END SUBROUTINE CalculateBodyAverage


  ! Create a permutation vector that maps elements (or nodes) to a smaller set
  ! of elements (or nodes) assuming periodicity in rotational angle.
  !----------------------------------------------------------------------------------
  SUBROUTINE RotationalPeriodicSumPerm(Solver, Mesh, angle, Perm, SumPerm, ElemField, IsSymmetric )
    TYPE(Solver_t) :: Solver
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: angle
    INTEGER :: Perm(:), SumPerm(:)
    LOGICAL :: ElemField, IsSymmetric
    
    INTEGER :: i,j,k,n,m,hits,nsym
    REAL(KIND=dp) :: x0,y0,r0,x1,y1,r1,phi0,phi1,dphi,phieps,reps,smax,maxdphi,phimin
    TYPE(Element_t), POINTER :: Element0, Element1
    INTEGER, POINTER :: Inds0(:), Inds1(:)
        
    phieps = 1.0e-3*angle
    reps = 1.0e-3
    maxdphi = 0.0_dp
    hits = 0    
    nsym = 0
    SumPerm = 0
    
    IF( ElemField ) THEN
      phimin = HUGE(phimin)
      n = 0

      DO j=1,Mesh % NumberOfBulkElements          
        Element1 => Mesh % Elements(j)
        Inds1 => Element1 % NodeIndexes
        IF(ANY(Perm(Inds1)==0)) CYCLE

        n = n+1        
        DO i=1,Element1 % TYPE % NumberOfNodes
          x1 = Mesh % Nodes % x(Inds1(i))
          y1 = Mesh % Nodes % y(Inds1(i))
          r1 = SQRT(x1*x1+y1*y1)

          phi1 = 180.0_dp*ATAN2(y1,x1)/PI
          IF(phi1 > 90.0) phi1 = phi1 - 180.0_dp
          phimin = MIN(phimin,phi1)
        END DO
      END DO
      CALL Info('CreatePeriodicSumPerm','Number of element in rotational piece: '//I2S(n),Level=15)
      WRITE(Message,'(A,ES12.5)') 'Offset of rotational piece: ',phimin
      CALL Info('CreatePeriodicSumPerm',Message,Level=10)
      
      DO i=1,Mesh % NumberOfBulkElements        
        Element0 => Mesh % Elements(i)
        Inds0 => Element0 % NodeIndexes

        n = Element0 % TYPE % NumberOfNodes
        x0 = SUM(Mesh % Nodes % x(Inds0)) / n
        y0 = SUM(Mesh % Nodes % y(Inds0)) / n
        r0 = SQRT(x0*x0+y0*y0)
        phi0 = 180.0_dp*ATAN2(y0,x0)/PI - phimin

        IF( IsSymmetric ) THEN
          phi0 = MODULO(phi0,2*angle)
        ELSE
          phi0 = MODULO(phi0,angle)
        END IF
        
        smax = MAXVAL(Mesh % Nodes % x(Inds0)) - MINVAL(Mesh % Nodes % x(Inds0)) &
            + MAXVAL(Mesh % Nodes % y(Inds0)) - MINVAL(Mesh % Nodes % y(Inds0))
        reps = 1.0e-3 * smax
        
        DO j=1,Mesh % NumberOfBulkElements          
          Element1 => Mesh % Elements(j)
          Inds1 => Element1 % NodeIndexes

          IF(Element1 % TYPE % NumberOfNodes /= n) CYCLE          
          IF(ANY(Perm(Inds1)==0)) CYCLE

          x1 = SUM(Mesh % Nodes % x(Inds1)) / n
          y1 = SUM(Mesh % Nodes % y(Inds1)) / n
          r1 = SQRT(x1*x1+y1*y1)
          IF(ABS(r1-r0) > reps ) CYCLE

          phi1 = 180.0_dp*ATAN2(y1,x1)/PI - phimin

          IF( IsSymmetric ) THEN
            phi1 = MODULO(phi1,2*angle)

            ! Periodic 2*angle ? 
            dphi = phi0-phi1                       
            IF(ABS(dphi) < phieps ) THEN
              SumPerm(i) = j            
              hits = hits+1
              EXIT
            END IF

            ! Test for symmetric hit
            dphi = 2*angle - (phi0+phi1)
            IF(ABS(dphi) < phieps ) THEN
              SumPerm(i) = -j            
              nsym = nsym+1
              hits = hits+1
              EXIT
            END IF
          ELSE
            phi1 = MODULO(phi1,angle)
            dphi = phi0-phi1            
            IF(ABS(dphi) < phieps ) THEN
              SumPerm(i) = j            
              hits = hits+1
              EXIT
            END IF
          END IF
            
          maxdphi = MAX(maxdphi,dphi)
        END DO
      END DO

      m = COUNT(SumPerm==0)
      CALL Info('CreatePeriodicSumPerm','Number of misses in rotational piece: '//I2S(m),Level=15)
      CALL Info('CreatePeriodicSumPerm','Elemental periodic perm with '//I2S(hits)//' hits',Level=10)            
    ELSE      
      DO i=1,Mesh % NumberOfNodes
        k = Perm(i)
        IF(k==0) CYCLE
        x0 = Mesh % Nodes % x(i)
        y0 = Mesh % Nodes % y(i)
        r0 = SQRT(x0*x0+y0*y0)
        phi0 = 180.0_dp*ATAN2(y0,x0)/PI
        DO j=1,Mesh % NumberOfNodes
          x1 = Mesh % Nodes % x(j)
          y1 = Mesh % Nodes % y(j)
          r1 = SQRT(x1*x1+y1*y1)
          phi1 = 180*ATAN2(y1,x1)/PI
          IF(ABS(r1-r0) < reps ) THEN
            IF(ABS(MODULO(phi0-phi1,dphi)) < phieps ) THEN
              SumPerm(j) = k            
              hits = hits+1
            END IF
          END IF
        END DO
      END DO
      CALL Info('CreatePeriodicSumPerm','Generated periodic sum perm with '//I2S(hits)//' hits')            
    END IF
      
  END SUBROUTINE RotationalPeriodicSumPerm
    
  

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
    

    CALL Info('MinimalElementalSet','Creating discontinuous subset from DG field',Level=5)

    ! Calculate size of permutation vector
    ALLOCATE( NodeVisited( Mesh % NumberOfNodes ) )
    NodeVisited = 0

    NULLIFY( SetPerm ) 
    k = 0
    DO i=1,Mesh % NumberOfBulkElements         
      Element => Mesh % Elements(i)
      k = k + Element % TYPE % NumberOfNodes
    END DO
    CALL Info('MinimalElementalSet','Maximum number of dofs in DG: '//I2S(k),Level=12)
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

    CALL Info('MinimalElementalSet','Independent dofs in elemental field: '//I2S(l),Level=7)
    CALL Info('MinimalElementalSet','Redundant dofs in elemental field: '//I2S(NoElimNodes),Level=7)     

  END FUNCTION MinimalElementalSet


  !> Calculate the reduced DG field given the reduction permutation.
  !> The permutation must be predefined. This may be called repeatedly
  !> for different variables. Optionally one may take average, or 
  !> a plain sum over the shared nodes. 
  !-------------------------------------------------------------------
  SUBROUTINE ReduceElementalVar( Mesh, Var, SetPerm, TakeAverage )

    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), TARGET :: Mesh
    INTEGER, TARGET :: SetPerm(:)
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
      CALL Info('ReduceElementalVar','Calculating reduced set average for: '&
          //TRIM(Var % Name), Level=7)
    ELSE
      CALL Info('ReduceElementalVar','Calculating reduced set sum for: '&
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
        
      m = SUM( SetCount ) 
      IF( m == 0 ) RETURN

      IF( TakeAverage ) THEN
        WHERE( SetCount > 0 ) SetSum = SetSum / SetCount
      END IF

      IF( dof == 1 ) THEN
        AveHits = 1.0_dp * SUM( SetCount ) / COUNT( SetCount > 0 )
        WRITE(Message,'(A,ES15.4)') 'Average number of hits: ',AveHits
        CALL Info('ReduceElementalVar',Message,Level=10)
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
    TYPE(Mesh_t), TARGET :: Mesh
    INTEGER, TARGET :: SetPerm(:)
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
        CALL Info('LumpedElementalVar','Lumped sum for component: '//I2S(dof),Level=6)
      END IF
      DO i=1,CurrentModel % NumberOfBodies
        WRITE(Message,'(A,ES15.4)') 'Body '//I2S(i)//' sum:',BodySum(i)
        CALL Info('LumpedElementalVar',Message,Level=10)
      END DO

    END DO

    DEALLOCATE( NodeVisited, BodySum )

  END SUBROUTINE LumpedElementalVar



!------------------------------------------------------------------------------
  SUBROUTINE SaveParallelInfo( Solver )
!------------------------------------------------------------------------------
   TYPE( Solver_t ), POINTER  :: Solver
!------------------------------------------------------------------------------    
   TYPE(ParallelInfo_t), POINTER :: ParInfo=>NULL()
   TYPE(ValueList_t), POINTER :: Params
   INTEGER :: i,j,k,n,maxnei
   LOGICAL :: Found, MeshMode, MatrixMode
   CHARACTER(*), PARAMETER :: Caller = "SaveParallelInfo"
   TYPE(Nodes_t), POINTER :: Nodes
   CHARACTER(:), ALLOCATABLE :: dumpfile
   
   Params => Solver % Values 

   MeshMode = ListGetLogical( Params,'Save Parallel Matrix Info',Found ) 
   MatrixMode = ListGetLogical( Params,'Save Parallel Mesh Info',Found ) 

   IF( .NOT. ( MeshMode .OR. MatrixMode ) ) RETURN

10 IF( MeshMode ) THEN
     CALL Info(Caller,'Saving parallel mesh info',Level=8 ) 
   ELSE
     CALL Info(Caller,'Saving parallel matrix info',Level=8 ) 
   END IF

   IF( MeshMode ) THEN
     ParInfo => Solver % Mesh % ParallelInfo
     Nodes => Solver % Mesh % Nodes
     dumpfile = 'parinfo_mesh.dat'
   ELSE
     ParInfo => Solver % Matrix % ParallelInfo
     dumpfile = 'parinfo_mat.dat'      
   END IF

   IF( .NOT. ASSOCIATED( ParInfo ) ) THEN
     CALL Warn(Caller,'Parallel info not associated!')
     RETURN
   END IF

   n = SIZE( ParInfo % GlobalDOFs )
   IF( n <= 0 ) THEN
     CALL Warn(Caller,'Parallel info size is invalid!')
     RETURN
   END IF

   ! memorize the maximum number of parallel neighbours
   maxnei = 0
   IF( ASSOCIATED( ParInfo % NeighbourList ) ) THEN
     DO i=1,n
       IF( ASSOCIATED( ParInfo % NeighbourList(i) % Neighbours ) ) THEN
         j = SIZE( ParInfo % NeighbourList(i) % Neighbours )
         maxnei = MAX( j, maxnei ) 
       END IF
     END DO
   END IF
   CALL Info(Caller,'Maximum number of parallel neighbours:'//I2S(maxnei))

   IF(ParEnv % PEs > 1) dumpfile = TRIM(dumpfile)//'.'//I2S(ParEnv % myPE)      
   CALL Info(Caller,'Saving parallel info to: '//TRIM(dumpfile),Level=8)

   OPEN(1,FILE=dumpfile, STATUS='Unknown')  
   DO i=1,n
     j = ParInfo % GlobalDOFs(i)
     IF( ParInfo % GInterface(i) ) THEN
       k = 1
     ELSE
       k = 0
     END IF
     WRITE(1,'(3I6)',ADVANCE='NO') i,j,k
     IF( ASSOCIATED( ParInfo % NeighbourList(i) % Neighbours ) ) THEN
       k = SIZE( ParInfo % NeighbourList(i) % Neighbours )
     ELSE
       k = 0
     END IF
     DO j=1,k
       WRITE(1,'(I6)',ADVANCE='NO')  ParInfo % NeighbourList(i) % Neighbours(j)
     END DO
     DO j=k+1,maxnei
       WRITE(1,'(I6)',ADVANCE='NO')  -1 
     END DO
     IF( MeshMode ) THEN
       WRITE(1,'(3ES12.3)',ADVANCE='NO') &
           Nodes % x(i), Nodes % y(i), Nodes % z(i)
     END IF
     WRITE(1,'(A)') ' ' ! finish the line
   END DO
   CLOSE(1)

   ! Redo with matrix if both modes are requested
   IF( MeshMode .AND. MatrixMode ) THEN
     MeshMode = .FALSE.
     GOTO 10
   END IF
   
   CALL Info(Caller,'Finished saving parallel info',Level=10)

!------------------------------------------------------------------------------
 END SUBROUTINE SaveParallelInfo
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetLagrangeIndexes( Mesh, LagN, Element, Indexes )  RESULT(L)
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    INTEGER :: LagN
    TYPE(Element_t), OPTIONAL, TARGET :: Element
    INTEGER, OPTIONAL :: Indexes(:)
    INTEGER :: L
!------------------------------------------------------------------------------    
    TYPE(Solver_t),  POINTER :: Solver
    TYPE(Element_t), POINTER :: Parent, Edge, Face
    LOGICAL :: OrientationsMatch
    LOGICAL :: EdgesActive, FacesActive
    LOGICAL :: Visited = .FALSE.

    INTEGER, PARAMETER :: MAX_LAGRANGE_NODES = 729

    INTEGER :: EdgeMap(2), FaceMap(4)
    INTEGER :: VTKTetraFaceMap(4,3)
    INTEGER :: VTKBrickFaceMap(6,4), BrickFaceOrdering(6)
    INTEGER :: Perm(MAX_LAGRANGE_NODES), TmpInd(MAX_LAGRANGE_NODES)
    INTEGER :: i,j,m,n0,e1,e2,f
    INTEGER :: nelem, nface, nedge, elemdim, thiselem
    INTEGER :: nelem_max, nface_max, nedge_max
    INTEGER :: ElemFamily, nsize
    INTEGER :: ElemType
    CHARACTER(*), PARAMETER :: Caller = 'GetLagrangeIndexes'

    SAVE Visited, nelem_max, nface_max, nedge_max, nsize, EdgesActive, &
        FacesActive, VTKTetraFaceMap, VTKBrickFaceMap, BrickFaceOrdering
!------------------------------------------------------------------------------
    
    IF (.NOT. Visited) THEN
      Visited = .TRUE.

      ! VTK's convention:
      VTKTetraFaceMap(1,:) = (/ 1,2,4 /)
      VTKTetraFaceMap(2,:) = (/ 3,4,2 /)
      VTKTetraFaceMap(3,:) = (/ 1,4,3 /)
      VTKTetraFaceMap(4,:) = (/ 1,3,2 /)

      VTKBrickFaceMap(1,:) = (/ 1,4,8,5 /)
      VTKBrickFaceMap(2,:) = (/ 2,3,7,6 /)
      VTKBrickFaceMap(3,:) = (/ 1,2,6,5 /)
      VTKBrickFaceMap(4,:) = (/ 4,3,7,8 /)
      VTKBrickFaceMap(5,:) = (/ 1,2,3,4 /)
      VTKBrickFaceMap(6,:) = (/ 5,6,7,8 /)
      BrickFaceOrdering = (/ 6,4,3,5,1,2 /)

      nedge_max = 0
      nface_max = 0
      nelem_max = 0

      DO i=1,Mesh % NumberOfBulkElements
        ElemFamily = Mesh % Elements(i) % TYPE % ElementCode / 100
        CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
        nedge_max = MAX(nedge, nedge_max)
        nface_max = MAX(nface, nface_max)
        nelem_max = MAX(nelem, nelem_max) 
      END DO
      
      EdgesActive = ASSOCIATED(Mesh % Edges)
      FacesActive = ASSOCIATED(Mesh % Faces)
      
      IF (.NOT. EdgesActive .AND. nedge_max > 0) CALL Warn(Caller, 'Mesh edges needed but not associated')
      IF (.NOT. FacesActive .AND. nface_max > 0) CALL Warn(Caller, 'Mesh faces needed but not associated')

      nsize = Mesh % NumberOfNodes + nelem_max * Mesh % NumberOfBulkElements + &
          nface_max * Mesh % NumberOfFaces + nedge_max * Mesh % NumberOfEdges

      nsize = nsize + Mesh % NumberOfBoundaryElements * MAX(nedge_max, nface_max)
    END IF

    ! If we don't have a specific element, then only return the total number which is sufficiently large
    ! in order to index all DOFs in the Lagrange mesh. 
    IF (.NOT. PRESENT(Element)) THEN
      l = nsize
      RETURN
    END IF
        
    ! The count of corner nodes:
    l = Element % TYPE % ElementCode / 100 
    IF( l >= 5 .AND. l <= 7 ) l = l-1             

    IF (PRESENT(Indexes)) THEN
      Indexes = 0
      Indexes(1:l) = Element % NodeIndexes(1:l)
    END IF
    ! Offset
    n0 = Mesh % NumberOfNodes

    IF(l>4) THEN
      ElemDim = 3
    ELSE IF(l>2) THEN
      ElemDim = 2
    ELSE
      ElemDim = 1
    END IF

    
    ! Number the additional edge nodes
    IF (EdgesActive ) THEN
      ElemFamily = Element % TYPE % ElementCode / 100
      CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)

      ! If this is a boundary element, we need to number it just as it would if it were an edge
      ! of a bulk element. 
      IF( ElemDim == 1 .AND. ASSOCIATED(Element % BoundaryInfo) ) THEN
        thiselem = 1        
        nedge = nelem
      ELSE
        thiselem = 0
      END IF
      
      DO i=1,MAX(thiselem,Element % TYPE % NumberOfEdges)
        IF(thiselem==1) THEN
          ! We use sneaky definitions here to be able to use rest of the edge indexing code.
          ! We want to use the edge indexing that has been generated for the edges of
          ! the parent element.
          Parent => Element % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED( Parent ) ) THEN
            Parent => Element % BoundaryInfo % Right
          END IF
          IF (.NOT. ASSOCIATED(Parent)) RETURN
          Edge => Find_Edge(Mesh,Parent,Element)
          EdgeMap = [1,2]
        ELSE
          f = i
          SELECT CASE(ElemFamily)
          CASE(2)
            CALL Error(Caller, '2D element is supposed to have elemental DOFs')
          CASE(3)
            EdgeMap = GetTriangleEdgeMap(i)
          CASE(4)
            EdgeMap = GetQuadEdgeMap(i)
          CASE(5)
            EdgeMap = GetTetraEdgeMap(i)
            IF (i == 3) THEN
              e1 = EdgeMap(2)
              e2 = EdgeMap(1)
              EdgeMap(1) = e1
              EdgeMap(2) = e2
            END IF
          CASE(6)
            EdgeMap = GetPyramidEdgeMap(i)
          CASE(7)
            EdgeMap = GetWedgeEdgeMap(i)
          CASE(8)
            ! It seems that VTK cell types 72 and 12/29 are not interchangeable:
            IF (LagN > 2) THEN
              ! The following is needed for 72:
              SELECT CASE(i)
              CASE(11)
                f = 12
              CASE(12)
                f = 11
              CASE DEFAULT
                CONTINUE
              END SELECT
            END IF
            EdgeMap = GetBrickEdgeMap(f)
          END SELECT
          Edge => Mesh % Edges(Element % EdgeIndexes(f))     
        END IF
        
        e1 = Edge % NodeIndexes(1)
        e2 = Edge % NodeIndexes(2)

        IF (e2 < e1) THEN
          OrientationsMatch = e1 == Element % NodeIndexes(EdgeMap(2))
        ELSE
          OrientationsMatch = e1 == Element % NodeIndexes(EdgeMap(1))
        END IF        
        
        ! Ensure the edge DOFs are listed in the right order:
        IF (OrientationsMatch) THEN
          DO j=1,nedge
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = n0 + nedge_max*(Edge % ElementIndex-1)+j
          END DO
        ELSE
          DO j=nedge,1,-1
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = n0 + nedge_max*(Edge % ElementIndex-1)+j
          END DO
        END IF
      END DO

      ! Nothing to be done here. This was boundary element that was exhausted.
      IF(thiselem==1) RETURN
      
      n0 = n0 + Mesh % NumberOfEdges * nedge_max      
    END IF

    ! Then number the additional face nodes
    IF (FacesActive) THEN
      
      SELECT CASE(Element % TYPE % ElementCode / 100)
      CASE(3,4)
        ! For 2D element only save the face if it is a boundary!
        IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
          Parent => Element % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED( Parent ) ) THEN
            Parent => Element % BoundaryInfo % Right
          END IF
          IF (.NOT. ASSOCIATED(Parent)) RETURN
          Face => Find_Face(Mesh,Parent,Element)
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          
          IF (nelem < 1) RETURN

          IF (ElemFamily == 4) THEN
            Perm = LagrangeQuadFacePermutation(Element % NodeIndexes(1:4), LagN)
          ELSE
            Perm(1:3) = LagrangeTriFacePermutation(Element % NodeIndexes(1:3), LagN)
          END IF

          IF (PRESENT(Indexes)) THEN
            DO j=1,nelem
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF
          ! Permute to create the final list of indices:
          DO j=1,nelem
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END IF
        RETURN          

      CASE(5)
        DO i=1,Element % Type % NumberOfFaces
          !
          ! Elmer has created its face indices by using face maps different from
          ! VTK's convention. Set f so that we can assign the right global indices
          ! to the face i according to VTK's convention.
          !
          IF (i == 4) THEN
            f = 1
          ELSE
            f = i+1
          END IF
          
          Face => Mesh % Faces(Element % FaceIndexes(f))          
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          ! test:
          !m = 0
          !DO j=1,3
          !  DO k=1,3
          !    IF (Face % NodeIndexes(j) == Element % NodeIndexes(VTKTetraFaceMap(i,k))) THEN
          !      m = m + 1
          !      EXIT
          !    END IF
          !  END DO
          !END DO
          !IF (m /= 3) CALL Fatal(Caller, 'Face is not identified correctly')

          Perm(1:3) = LagrangeTriFacePermutation(Element % NodeIndexes(VTKTetraFaceMap(i,1:3)), LagN)
          
          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF

          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO

      CASE(6)
        ! The quad face:
        Face => Mesh % Faces(Element % FaceIndexes(1))
        CALL LagrangeDOFCount(4, LagN, nedge, nface, nelem)
        nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D
        FaceMap = GetPyramidFaceMap(1)

        IF (nface > 1) THEN
          CALL Fatal(Caller, 'For pyramids Lagrange Element Degree < 3 supported currently')
        END IF

        DO j=1,nface
          l = l + 1
          IF (PRESENT(Indexes)) Indexes(l) = n0 + nface_max*(Face % ElementIndex-1) + j
        END DO

        ! TO DO: Index triangular faces for degrees p > 3

      CASE(7)
        ! Triangular faces:
        DO f=1,2
          Face => Mesh % Faces(Element % FaceIndexes(f))
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          IF (nface < 1) CYCLE

          FaceMap = GetWedgeFaceMap(f)
          Perm(1:3) = LagrangeTriFacePermutation(Element % NodeIndexes(FaceMap(1:3)), LagN)
          
          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF

          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO

        ! Quad faces:
        DO f=3,5
          Face => Mesh % Faces(Element % FaceIndexes(f))          
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          IF (nface < 1) CYCLE

          FaceMap = GetWedgeFaceMap(f)
          Perm = LagrangeQuadFacePermutation(Element % NodeIndexes(FaceMap(1:4)), LagN)

          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF

          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO

      CASE(8)
        DO i=1,Element % Type % NumberOfFaces 
          f = BrickFaceOrdering(i)

          Face => Mesh % Faces(Element % FaceIndexes(f))          
          ElemFamily = Face % TYPE % ElementCode / 100
          CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)
          nface = nelem ! The number of elementwise DOFs in 2D gives the count of face DOFs in 3D

          IF (nface < 1) CYCLE

          Perm = LagrangeQuadFacePermutation(Element % NodeIndexes(VTKBrickFaceMap(i,1:4)), LagN)

          IF (PRESENT(Indexes)) THEN
            DO j=1,nface
              TmpInd(j) = n0 + nface_max*(Face % ElementIndex-1) + j
            END DO
          END IF
          ! Permute to create the final list of indices:
          DO j=1,nface
            l = l + 1
            IF (PRESENT(Indexes)) Indexes(l) = TmpInd(Perm(j))
          END DO
        END DO
      END SELECT
      
      n0 = n0 + Mesh % NumberOfFaces * nface_max
    END IF

    ! Then number the additional internal nodes (never shared)
    ElemFamily = Element % TYPE % ElementCode / 100
    CALL LagrangeDOFCount(ElemFamily, LagN, nedge, nface, nelem)    
    DO j=1,nelem
      l = l + 1
      IF (PRESENT(Indexes)) Indexes(l) = n0 + nelem_max*(Element % ElementIndex-1) + j
    END DO

  CONTAINS
    ! 
    ! A subroutine for returning the maximal number of interior nodes associated with
    ! the element edges, faces and the volume in the Lagrange interpolation of degree p
    !
    SUBROUTINE LagrangeDOFCount(Family, p, nedge, nface, nelem)
      INTEGER, INTENT(IN) :: Family, p
      INTEGER, INTENT(OUT) :: nedge, nface, nelem
      
      INTEGER :: m

      m = p - 1
      nelem = 0
      nface = 0
      nedge = 0

      IF (Family == 1) RETURN
      
      SELECT CASE(Family)
      CASE(2)
        nelem = m
      CASE(3)
        nelem = m*(m-1)/2
        nedge = m
      CASE(4)
        nelem = m*m
        nedge = m
      CASE(5)
        nelem = m*(m-1)*(m-2)/6
        nface = m*(m-1)/2
        nedge = m
      CASE(6)
        nedge = m
        nface = m*m ! the maximum is determined by quad faces
        IF (p > 1) THEN
          IF (p==2) THEN
            nelem = 1
          ELSE
            CALL Fatal('LagrangeDOFCount', 'Cannot handle pyramids of degree > 2')
          END IF
        END IF
      CASE(7)
        nedge = m
        nface = m*m ! the maximum is determined by quad faces
        nelem = m*(m-1)/2*m
      CASE(8)
        nelem = m*m*m
        nface = m*m
        nedge = m
      CASE DEFAULT          
        CALL Fatal('LagrangeDOFCount', 'Unknown element family') 
      END SELECT
    END SUBROUTINE LagrangeDOFCount

    !
    ! A function to generate a permutation vector for indexing nodes on quad faces
    !
    FUNCTION LagrangeQuadFacePermutation(FaceNodes, p) RESULT(Perm)
      INTEGER, INTENT(IN) :: FaceNodes(4)
      INTEGER, INTENT(IN) :: p       ! the order of Lagrange interpolation
      INTEGER, PARAMETER :: MAX_LAGRANGE_NODES = 729
      INTEGER :: Perm(MAX_LAGRANGE_NODES)
      INTEGER :: AllIndices((p-1)**2)
      INTEGER :: i, j, n, i0, MinEntryInd(1)

      SELECT CASE(p)
      CASE(2)
        Perm = 0
        Perm(1) = 1

      CASE DEFAULT
        !
        ! We have 4 x 2 permutation patterns. Create a permutation
        ! vector to alter the default ordering in each case. The first face
        ! index is assigned to the node which is closest to the face corner A
        ! having the smallest global index. The next indices are created in 
        ! the direction of the face edge AB, with B the smallest possible 
        ! global index.
        !
        Perm = 0
        n = (p-1)**2
        DO i=1,n
          AllIndices(i) = i
        END DO

        MinEntryInd = MINLOC(FaceNodes(1:4))
        SELECT CASE(MinEntryInd(1))
        CASE(1)
          IF (FaceNodes(4) < FaceNodes(2)) THEN
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              Perm(i0+1:i0+p-1) = AllIndices(i:n:p-1)
            END DO
          ELSE
            Perm(1:n) = AllIndices(1:n)
          END IF

        CASE(2)
          IF (FaceNodes(3) < FaceNodes(1)) THEN
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(p-i+(j-1)*(p-1))
              END DO
            END DO
          ELSE
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(i0+p-j)
              END DO
            END DO
          END IF          

        CASE(3)
          IF (FaceNodes(4) < FaceNodes(2)) THEN
            DO i=1,n
              Perm(i) = AllIndices(n+1-i)
            END DO
          ELSE
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(n+1-i-(j-1)*(p-1))
              END DO
            END DO
          END IF

        CASE(4)
          IF (FaceNodes(1) < FaceNodes(3)) THEN
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(n-p+1+i-(j-1)*(p-1))
              END DO
            END DO
          ELSE
            DO i=1,p-1
              i0 = (i-1)*(p-1)
              DO j=1,p-1
                Perm(i0+j) = AllIndices(n-i*(p-1)+j)
              END DO
            END DO
          END IF
        END SELECT

      END SELECT
    END FUNCTION LagrangeQuadFacePermutation

    !
    ! A function to generate a permutation vector for indexing nodes on triangular faces
    !
    FUNCTION LagrangeTriFacePermutation(FaceNodes, p) RESULT(Perm)
      INTEGER, INTENT(IN) :: FaceNodes(3)
      INTEGER, INTENT(IN) :: p       ! the order of Lagrange interpolation
      INTEGER :: Perm(3)

      INTEGER :: MinEntryInd(1)

      SELECT CASE(p)
      CASE(3)
        Perm = 0
        Perm(1) = 1

      CASE(4)
        !
        ! We have 3 x 2 permutation patterns. Create a permutation
        ! vector to alter the default ordering in each case. The first face
        ! index is assigned to the node which is closest to the face corner A
        ! having the smallest global index. The next indices are created in 
        ! the direction of the face edge AB, with B the smallest possible 
        ! global index.
        !
        Perm = 0

        MinEntryInd = MINLOC(FaceNodes(1:3))
        SELECT CASE(MinEntryInd(1))
        CASE(1)
          IF (FaceNodes(3) < FaceNodes(2)) THEN
            Perm = (/ 1,3,2 /)
          ELSE
            Perm = (/ 1,2,3 /)
          END IF

        CASE(2)
          IF (FaceNodes(3) < FaceNodes(1)) THEN
            Perm = (/ 2,3,1 /)
          ELSE
            Perm = (/ 2,1,3 /)
          END IF          

        CASE(3)
          IF (FaceNodes(1) < FaceNodes(2)) THEN
            Perm = (/ 3,1,2 /)
          ELSE
            Perm = (/ 3,2,1 /)
          END IF
        END SELECT

      CASE DEFAULT
        CALL Fatal('LagrangeTriFacePermutation', &
            'For triangular faces Lagrange Element Degree < 5 supported currently')

      END SELECT
    END FUNCTION LagrangeTriFacePermutation



!------------------------------------------------------------------------------
  END FUNCTION GetLagrangeIndexes
!------------------------------------------------------------------------------   


 !> Find a representative DG index for a node index. Note that
 !> there may be several possibilities and this is just one of them.
 !------------------------------------------------------------------  
   FUNCTION NodeToDGIndex(Mesh,nodeind) RESULT ( dgind )

    TYPE(Mesh_t) :: Mesh
    INTEGER :: nodeind
    INTEGER :: dgind
    
    INTEGER :: i,j,t
    TYPE(Element_t), POINTER :: Element

    dgind = 0

    IF(nodeind < 1 ) THEN
      CALL Warn('NodeToDGIndex','Cannot find DG index for too small node index!')
      RETURN
    END IF
    IF(nodeind > Mesh % NumberOfNodes ) THEN
      CALL Warn('NodeToDGIndex','Cannot find DG index for too large node index!')
      RETURN
    END IF         
    
    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      DO i = 1,Element % TYPE % NumberOfNodes          
        IF( Element % NodeIndexes(i) == nodeind ) THEN
          IF(.NOT. ASSOCIATED( Element % DGIndexes ) ) THEN
            CALL Fatal('NodeToDGIndex','There are no DG indexes!')
          END IF
          dgind = Element % DGIndexes(i)
          EXIT
        END IF
      END DO
      IF(dgind > 0 ) EXIT
    END DO
    
  END FUNCTION NodeToDGIndex



!------------------------------------------------------------------------------
!> Split a mesh at zero levelset by adding new nodes at the interface.
!> The idea is to be able to better represent shapes that are not initially
!> presented by body fitted finite element mesh. 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------


  !> Create interface boundaries consisting of edges defined by the intersection of two higher
  !> dimensional boundaries. This may be useful for 3D meshes where 1D meshes have not been
  !> create in advance.
  !-------------------------------------------------------------------



  !> Sometimes the mesh includes boundaries but it is annoyingly time-consuming to tag
  !> them by hand as their numbers are not known. Then an alternative is to use some simple
  !> rules to detect the existing boundaries and tag them with to the boundary that fulfills
  !> the detection rule.
  !-------------------------------------------------------------------


  !> Sometimes the mesh includes bodies different from how we would like to resolve
  !> the equations. Then an alternative is to use some simple rule to redefine the
  !> body index.
  !-------------------------------------------------------------------

  ! This is separated from the general UnitSegmentDivision since it can be used
  ! in other places as well. Note that w should range from 0 to n.
  !----------------------------------------------------------------------------
  SUBROUTINE GeometricUnitDivision(w, n, q)
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    INTEGER :: n
    REAL(KIND=dp) :: q

    INTEGER :: i,j
    REAL(KIND=dp) :: r,h1

    IF( n < 1 ) THEN
      CALL Fatal('GeometricUnitDivision','Cannot create division for '//I2S(n)//' element!')
    ELSE IF( ( ABS(ABS(q)-1.0_dp) < 1.0d-6 ) .OR. (q < 0.0_dp .AND. n <= 2) .OR. n==1) THEN
      CALL Info('GeometricUnitDivision','Creating linear division',Level=8)
      DO i=0,n
        w(i) = i/(1._dp * n)
      END DO
    ELSE
      CALL Info('GeometricUnitDivision','Creating geometric division',Level=8)
      IF( q > 0.0_dp ) THEN
        r = q**(1.0_dp/(n-1))
        h1 = (1-r)/(1-r**n)
        w(0) = 0.0_dp
        DO i=1,n-1
          w(i) = h1 * (1-r**i)/(1-r)
        END DO
        w(n) = 1.0_dp
      ELSE
        q = -q
        IF(MODULO(n,2) == 0) THEN
          r = q**(1.0_dp/(n/2-1))
          h1 = 0.5_dp*(1-r)/(1-r**(n/2))
        ELSE
          r = q**(1.0_dp/((n-1)/2))
          h1 = 0.5_dp / ( (1-r**((n+1)/2))/(1-r) - 0.5_dp * r**((n-1)/2))
        END IF

        w(0) = 0.0_dp
        DO i=1,n
          IF( i <= n/2 ) THEN
            w(i) = h1 * (1-r**i)/(1-r)
          ELSE
            w(i) = 1.0_dp -  h1 * (1-r**(n-i))/(1-r)
          END IF
        END DO
        w(n) = 1.0_dp
      END IF
    END IF

  END SUBROUTINE GeometricUnitDivision


  ! This is separated from the general UnitSegmentDivision since it can be used
  ! in other places as well. Note that w should range from 0 to n.
  !----------------------------------------------------------------------------
  SUBROUTINE FunctionUnitDivision(w, n, FunName, FunList )
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    INTEGER :: n
    CHARACTER(:), ALLOCATABLE :: FunName
    TYPE(ValueList_t), POINTER :: FunList

    INTEGER :: i,j,iter,maxiter
    REAL(KIND=dp) :: r,h1,hn,minhn,err_eps,err,xn
    REAL(KIND=dp), ALLOCATABLE :: wold(:),h(:)

    CALL Info('FunctionUnitDivision','Creating functional division: '//TRIM(FunName),Level=5)

    IF( n < 1 ) THEN
      CALL Fatal('GeometricUnitDivision','Cannot create division for '//I2S(n)//' element!')
    END IF

    ! Initial guess is an even distribution
    DO i=0,n
      w(i) = i/(1._dp * n)
    END DO
    IF(n == 1 ) RETURN

    ALLOCATE( wold(0:n),h(1:n))
    wold = w

    ! parameters that determine the accuracy of the iteration
    maxiter = 10000
    err_eps = 1.0d-6

    ! Iterate to have a density distribution
    !---------------------------------------
    DO iter=1,maxiter

      minhn = HUGE(minhn)
      wold = w

      DO i=1,n
        xn = (w(i)+w(i-1))/2.0_dp
        minhn = MIN( minhn, w(i)-w(i-1) )
        h(i) = ListGetFun( FunList, FunName, xn )
        IF( h(i) < EPSILON( h(i) ) ) THEN
          CALL Fatal('FunctionUnitDivision','Given value for h(i) was negative!')
        END IF
      END DO

      DO i=1,n-1
        w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
      END DO
      DO i=n-1,1,-1
        w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
      END DO

      err = MAXVAL( ABS(w-wold))/minhn

      IF( err < err_eps ) THEN
        WRITE( Message, '(A,I0,A)') 'Convergence obtained in ',iter,' iterations'
        CALL Info('FunctionUnitDivision', Message, Level=9 )
        EXIT
      END IF
    END DO

    IF( iter > maxiter ) THEN
      CALL Warn('FunctionUnitDivision','No convergence obtained for the unit mesh division!')
    END IF
  END SUBROUTINE FunctionUnitDivision


!------------------------------------------------------------------------------
!> Create node distribution for a unit segment x \in [0,1] with n elements
!> i.e. n+1 nodes. There are different options for the type of distribution.
!> 1) Even distribution
!> 2) Geometric distribution
!> 3) Arbitrary distribution determined by a functional dependence
!> Note that the 3rd algorithm involves iterative solution of the nodal
!> positions and is therefore not bullet-proof.
!------------------------------------------------------------------------------
  SUBROUTINE UnitSegmentDivision( w, n, ExtList )
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    INTEGER :: n
    TYPE(ValueList_t), POINTER, OPTIONAL :: ExtList
    !---------------------------------------------------------------
    INTEGER :: i
    REAL(KIND=dp) :: q
    CHARACTER(:), ALLOCATABLE :: FunName
    LOGICAL :: Found, GotRatio, GotFun
    TYPE(ValueList_t), POINTER :: ParList

    IF( PRESENT( ExtList ) ) THEN
      ParList => ExtList
    ELSE
      ParList => CurrentModel % Simulation
    END IF

    DO i=1,2
      IF(i==1) THEN
        FunName = 'Extruded Mesh Density'
      ELSE
        FunName = '1D Mesh Density'
      END IF
      GotFun = ListCheckPresent( ParList, FunName )
      IF(GotFun) EXIT
    END DO

    IF( GotFun ) THEN
      CALL FunctionUnitDivision(w,n,FunName,ParList)
    ELSE
      q = ListGetConstReal( ParList,'Extruded Mesh Ratio',GotRatio)
      IF(.NOT. GotRatio) q = ListGetConstReal( ParList,'1D Mesh Ratio',GotRatio)
      IF(.NOT. GotRatio) q = 1.0_dp
      CALL GeometricUnitDivision(w,n,q)
    END IF

    IF(InfoActive(9)) THEN
      CALL Info('UnitSegmentDivision','Mesh division ready' )
      DO i=0,n
        WRITE( Message, '(A,I0,A,ES12.4)') 'w(',i,') : ',w(i)
        CALL Info('UnitSegmentDivision', Message )
      END DO
    END IF

  END SUBROUTINE UnitSegmentDivision

!------------------------------------------------------------------------------
END MODULE MeshBasics
!------------------------------------------------------------------------------

!> \}

