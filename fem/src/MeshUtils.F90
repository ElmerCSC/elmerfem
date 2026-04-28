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

MODULE MeshUtils

    USE ElementDescription
    USE BandwidthOptimize
    USE Interpolation
    USE ParallelUtils
    USE Lists
    USe ListMatrix
    USE MeshAllocations
    USE ElementUtils, ONLY : mGetBoundaryIndexesFromParent, &
        NormalDirection, CreateMatrix, TangentDirections, &
        FreeMatrix, Find_Face, Find_Edge
    USE MortarUtils, ONLY : MarkHaloNodes, GeneratePeriodicProjectors
    USE GeometryFitting, ONLY : CylinderFit, SphereFit, TorusFit
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

     LOGICAL :: Found
     TYPE(Element_t) :: Element

     ! Sanity check to avoid memory leaks
     IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
        ALLOCATE(Element % PDefs, STAT=istat)
        IF ( istat /= 0) CALL Fatal('AllocatePDefinitions','Unable to allocate memory')
     ELSE
       CALL Info('AllocatePDefinitions','P element definitions already allocated',Level=32)
     END IF

     ! Initialize fields
     Element % PDefs % P = 0 
     Element % PDefs % TetraType = 0
     Element % PDefs % isEdge = .FALSE.
     Element % PDefs % localNumber = 0
     Element % PDefs % GaussPoints = 0

     Element % PDefs % Serendipity = ListGetLogical( CurrentModel % Simulation, &
           'Serendipity P elements', Found )
     IF(.NOT.Found) Element % PDefs % Serendipity = .TRUE.
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
     Element % BoundaryInfo % Constraint =  0
     Element % BoundaryInfo % RadiationFactors => NULL()

!------------------------------------------------------------------------------
   END SUBROUTINE AllocateBoundaryInfo
!------------------------------------------------------------------------------




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
    TYPE(Mesh_t), POINTER :: Mesh
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
     INTEGER :: i,j,t,n,t1,t2,right,hits,nact,npass,nfalse
     TYPE(Element_t), POINTER :: Parent, Element
     
     t1 = Mesh % NumberOfBulkElements
     t2 = Mesh % NumberOfBoundaryElements
     nfalse = 0
     
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
         CALL Warn('DropFalseParents','Boundary element '//I2S(t)//' no longer has parents with same indexes!')
       END IF

       nfalse = nfalse + npass
     END DO

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
   INTEGER, POINTER :: DiscontPerm(:)

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



 !> Fortran reader for Elmer ascii and binary mesh file format.
 !> The ascii format is tried out first, if not success, binary is followed. 
 !> This is a Fortran replacement for the old C++ eio library. 
 !------------------------------------------------------------------------
 SUBROUTINE ElmerMeshReader(Step, PMesh, MeshNamePar, ThisPe, NumPEs, IsParallel )

   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel

   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: PrevStep=0, iostat
   INTEGER, PARAMETER :: FileUnit = 10
   INTEGER :: i,j,k,n,BaseNameLen, SharedNodes = 0, mype = 0, numprocs = 0
   INTEGER, POINTER :: NodeTags(:), ElementTags(:), LocalPerm(:)
   INTEGER :: MinNodeTag = 0, MaxNodeTag = 0, istat
   LOGICAL :: ElementPermutation=.FALSE., NodePermutation=.FALSE., Parallel, &
       PseudoParallel, Found
   CHARACTER(:), ALLOCATABLE :: BaseName, FileName


   SAVE PrevStep, BaseName, BaseNameLen, Mesh, mype, Parallel, &
       NodeTags, ElementTags, LocalPerm, PseudoParallel

   CALL Info('ElmerMeshReader','Performing step: '//I2S(Step),Level=8)

   IF( Step - PrevStep /= 1 ) THEN
     CALL Fatal('ElmerMeshReader','The routine should be called in sequence: '// &
         I2S(PrevStep)//' : '//I2S(Step) )
   END IF
   PrevStep = Step
   IF( PrevStep == 6 ) PrevStep = 0 

   IF( Step == 1 ) THEN
     IF(.NOT. PRESENT( MeshNamePar ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give MeshNamePar!')
     END IF
     BaseName = TRIM( MeshNamePar ) 
     IF(.NOT. PRESENT( PMesh ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give PMesh!')
     END IF
     Mesh => PMesh
     IF(.NOT. PRESENT( ThisPe ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give ThisPe!')
     END IF
     mype = ThisPe 
     IF(.NOT. PRESENT( NumPEs) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give NumPEs!')
     END IF
     numprocs = NumPEs
     IF(.NOT. PRESENT( IsParallel ) ) THEN
       CALL Fatal('ElmerMeshReader','When calling in mode one give IsParallel!')
     END IF
     Parallel = IsParallel

     PseudoParallel = .FALSE.
     IF(.NOT. Parallel ) THEN
       IF( ParEnv % PEs > 1 ) THEN
         PseudoParallel = ListGetLogical(CurrentModel % Simulation,'Enforce Parallel',Found ) 
         IF(.NOT. Found ) PseudoParallel = ListGetLogicalAnySolver(CurrentModel,'Enforce Parallel')
       END IF
     END IF
     
     i = LEN_TRIM(MeshNamePar)
     DO WHILE(MeshNamePar(i:i) == CHAR(0))
       i=i-1
     END DO
     BaseNameLen = i
     CALL Info('ElmerMeshReader','Base mesh name: '//TRIM(MeshNamePar(1:BaseNameLen)))
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
     IF( PseudoParallel ) THEN
       CALL InitPseudoParallel()
     ELSE
       CALL InitParallelInfo()
       CALL ReadSharedFile()
     END IF
       
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
          '/partitioning.'//I2S(numprocs)//&
           '/part.'//I2S(mype+1)//'.header'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.header'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadHeaderFile','Reading header info from file: '//TRIM(FileName),Level=10)
     END IF

     READ(FileUnit,*,IOSTAT=iostat) Mesh % NumberOfNodes, &
         Mesh % NumberOfBulkElements,&
         Mesh % NumberOfBoundaryElements
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not read header 1st line in file: '//TRIM(FileName))
     END IF

     Types = 0
     CountByType = 0
     READ(FileUnit,*,IOSTAT=iostat) TypeCount
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadHeaderFile','Could not read the type count in file: '//TRIM(FileName))
     END IF
     DO i=1,TypeCount
       READ(FileUnit,*,IOSTAT=iostat) Types(i),CountByType(i)
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadHeaderFile','Could not read type count '&
             //I2S(i)//'in file: '//TRIM(FileName))
       END IF
     END DO

     IF( Parallel ) THEN
       READ(FileUnit,*,IOSTAT=iostat) SharedNodes
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadHeaderFile','Could not read shared nodes in file: '//TRIM(FileName))
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

     !USE iso_c_binding
     REAL(c_double) :: Coords(3)
     REAL(c_float) :: SCoords(3)
     INTEGER :: NodeTag
     LOGICAL :: Binary, singlePrec

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//I2S(numprocs)//&
           '/part.'//I2S(mype+1)//'.nodes'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.nodes'
     END IF

     Binary = .FALSE.
     SinglePrec = .FALSE.
     
     OPEN( Unit=FileUnit, File=FileName, STATUS='old', ACTION='read', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       ! ascii file was not successfull, try with binary.
       Binary = .TRUE.
       OPEN( Unit=FileUnit, File=TRIM(FileName)//".bin", FORM='unformatted', &
           ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
       IF(iostat /= 0 ) THEN         
         SinglePrec = .TRUE.
         OPEN( Unit=FileUnit, File=TRIM(FileName)//".sbin", FORM='unformatted', &
             ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
       END IF
     END IF
     
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadNodesFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadNodesFile','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ALLOCATE( NodeTags(Mesh % NumberOfNodes ) ) 
     NodeTags = 0

     NodePermutation = .FALSE.
     DO j = 1, Mesh % NumberOfNodes
       IF(SinglePrec) THEN
         READ(FileUnit,IOSTAT=iostat) NodeTag, SCoords
         Coords = SCoords
       ELSE IF(Binary) THEN
         READ(FileUnit,IOSTAT=iostat) NodeTag, Coords
       ELSE
         READ(FileUnit,*,IOSTAT=iostat) NodeTag, k, Coords
       END IF
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadNodesFile','Problem load node '//I2S(j)//' in file: '//TRIM(Filename))
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
     LOGICAL :: halo, Binary 


     CALL AllocateVector( ElementTags, Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements, 'ReadElementsFile')   
     ElementTags = 0
     ElementPermutation = .FALSE.

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)// &
          '/partitioning.'//I2S(numprocs)//&
             '/part.'//I2S(mype+1)//'.elements'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.elements'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='old', iostat=IOSTAT )
     IF( iostat == 0 ) THEN
       Binary = .FALSE.
     ELSE
       ! ascii file was not successfull, try with binary.
       Binary = .TRUE.       
       OPEN( Unit=FileUnit, File=TRIM(FileName)//".bin", FORM='unformatted', &
           ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
     END IF

     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadElementsFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadElementsFile','Reading bulk elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=1,Mesh % NumberOfBulkElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadElementsFile','Element '//I2S(i)//' not associated!')
       END IF

       IF(Binary) THEN
         READ(FileUnit,IOSTAT=iostat) Tag, PartN, body, elemtype
       ELSE
         READ(FileUnit, '(a)', IOSTAT=iostat) str
         IF( iostat /= 0 ) THEN
           CALL Fatal('ReadElementsFile','Could not read start of element entry: '//I2S(j))
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
       END IF
         
       ElementTags(j) = tag
       IF( j /= tag ) ElementPermutation = .TRUE.             
       Element % ElementIndex = j
       Element % BodyId = body

       IF( partn > 0 ) THEN
         Element % PartIndex = partn-1
       ELSE
         Element % PartIndex = mype
       END IF

       Element % TYPE => GetElementType(ElemType)

       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadElementsFile','Element of type '&
             //I2S(ElemType)//' could not be associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       CALL AllocateVector( Element % NodeIndexes, n )

       IF( Binary ) THEN
         READ(FileUnit,IOSTAT=iostat) Element % NodeIndexes(1:n)
       ELSE
         IF( nread < n + ioffset + 3 ) THEN
           CALL Fatal('ReadElementsFile','Line '//I2S(j)//' does not contain enough entries')
         END IF
         Element % NodeIndexes(1:n) = IVals(4+ioffset:nread)
       END IF
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
     LOGICAL :: halo, Binary

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//I2S(numprocs)//&
           '/part.'//I2S(mype+1)//'.boundary'
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
       CALL AllocateVector( LocalEPerm, MaxEIndex - MinEIndex + 1, 'ReadBoundaryFile' )
       LocalEPerm = 0
       DO i=1,Mesh % NumberOfBulkElements
         LocalEPerm( ElementTags(i) - MinEIndex + 1 ) = i
       END DO
     ELSE
       MinEIndex = 1 
       MaxEIndex = Mesh % NumberOfBulkElements
     END IF


     OPEN( Unit=FileUnit, File=FileName, STATUS='old', iostat=IOSTAT )

     IF( iostat == 0 ) THEN
       Binary = .FALSE.
     ELSE
       ! ascii file was not successfull, try with binary.
       Binary = .TRUE.       
       OPEN( Unit=FileUnit, File=TRIM(FileName)//".bin", FORM='unformatted', &
           ACCESS = 'stream', STATUS='old', ACTION='read', IOSTAT = iostat )
     END IF

     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadBoundaryFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadBoundaryFile','Reading boundary elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadBoundaryFile','Element '//I2S(i)//' not associated!')
       END IF

       IF(Binary) THEN
         READ(FileUnit,IOSTAT=iostat) Tag, PartN, bndry, left, right, elemtype
       ELSE
         READ(FileUnit, '(a)', IOSTAT=iostat) str
         IF( iostat /= 0 ) THEN
           CALL Fatal('ReadBoundaryFile','Could not read boundary element entry: '//I2S(j))
         END IF
         nread = read_ints(str,ivals,halo)
         
         tag = ivals(1)
         ElementTags(j) = tag
         
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
       END IF
         
       Element % ElementIndex = j
       Element % TYPE => GetElementType(ElemType)
       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadBoundaryFile','Element of type '//I2S(ElemType)//'could not be associated!')
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

       IF( binary ) THEN
         READ(FileUnit,IOSTAT=iostat) Element % NodeIndexes(1:n)
       ELSE
         IF( nread < 5 + n + ioffset ) THEN
           CALL Fatal('ReadBoundaryFile','Line '//I2S(j)//' does not contain enough entries')
         END IF
         Element % NodeIndexes(1:n) = Ivals(6+ioffset:nread)
       END IF
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
       CALL Info('PermuteNodeNumbering','Performing node mapping',Level=6)

       MinNodeTag = MINVAL( NodeTags )
       MaxNodeTag = MAXVAL( NodeTags )

       CALL AllocateVector( LocalPerm, MaxNodeTag-MinNodeTag+1, 'PermuteNodeNumbering' )
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
       CALL Info('PermuteNodeNumbering','Node mapping is continuous',Level=8)
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


     ! This also for serial runs ...
     DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i)
     END DO

     IF(.NOT. Parallel ) RETURN

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
     IF (istat /= 0) CALL Fatal('InitParallelInfo', 'Unable to allocate NeighbourList array.')

     DO i=1,n
       NULLIFY( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % GInterface, n, 'InitParallelInfo')
     Mesh % ParallelInfo % GInterface = .FALSE.       

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
       '/partitioning.'//I2S(numprocs)//&
         '/part.'//I2S(mype+1)//'.shared'

     OPEN( Unit=FileUnit, File=FileName, STATUS='old', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadSharedFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('ReadSharedFile','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ! This loop could be made more effective, for example
     ! by reading tags and nparts to a temporal vector
     ! The operation using the str takes much more time.
     !-----------------------------------------------------
     DO i=1,SharedNodes          
       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadSharedFile','Could not read shared nodes entry: '//I2S(i))
       END IF
       nread = read_ints(str,ivals,halo)

       tag = ivals(1)
       npart = ivals(2)       

       k = LocalPerm( tag-MinNodeTag+1 )
       Mesh % ParallelInfo % GInterface(k) = .TRUE.
       CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(k) % Neighbours,npart)

       IF( nread < 2 + npart ) THEN
         CALL Fatal('ReadSharedFile','Line '//I2S(j)//' does not contain enough entries')
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


   ! Initialize parallel info for pseudo parallel meshes
   !-------------------------------------------------------
   SUBROUTINE InitPseudoParallel()

     INTEGER, POINTER :: TmpGlobalDofs(:)

     ! This also for serial runs ...
     n = ParEnv % MyPe * Mesh % NumberOfBulkElements

     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i) + n
     END DO

     n = Mesh % NumberOfNodes + &
         Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
         Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
         Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements

     ALLOCATE( TmpGlobalDOFs(n) )
     TmpGlobalDOFs = 0
     TmpGlobalDOFs(1:Mesh % NumberOfNodes) = &
         Mesh % ParallelInfo % GlobalDOFs(1:Mesh % NumberOfNodes) + n
     DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs ) 
     Mesh % ParallelInfo % GlobalDofs => TmpGlobalDofs
     
     ALLOCATE(Mesh % ParallelInfo % NeighbourList(n), STAT=istat)
     IF (istat /= 0) CALL Fatal('InitParallelInfo', 'Unable to allocate NeighbourList array.')
     
     DO i=1,n
       ALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) )
       Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % MyPe
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % GInterface, n, 'InitParallelInfo')
     Mesh % ParallelInfo % GInterface = .FALSE.       

   END SUBROUTINE InitPseudoParallel

   
 END SUBROUTINE ElmerMeshReader

 !> An interface over potential mesh loading strategies. 
 !----------------------------------------------------------------- 
 SUBROUTINE LoadMeshStep( Step, PMesh, MeshNamePar, ThisPe, NumPEs,IsParallel ) 
   
   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel

   ! Currently only one strategy to get the mesh is implemented 
   ! but there could be others.
   !
   ! This has not yet been tested in parallel and for sure
   ! it does not work for halo elements. 
   !-----------------------------------------------------------------
   CALL ElmerMeshReader( Step, PMesh, MeshNamePar, ThisPe, NumPEs, IsParallel ) 

 END SUBROUTINE LoadMeshStep

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

 ! Set the mesh dimension by studying the coordinate values.
 ! This could be less conservative also...
 !------------------------------------------------------------------------------    
 SUBROUTINE SetMeshDimension( Mesh )
   TYPE(Mesh_t), POINTER :: Mesh
   
   REAL(KIND=dp) :: x, y, z
   LOGICAL :: C(3)
   INTEGER :: i
   
   IF( Mesh % NumberOfNodes == 0 ) RETURN

   ! Compare value to some node, why not the 1st one
   x = Mesh % Nodes % x(1)
   y = Mesh % Nodes % y(1)
   z = Mesh % Nodes % z(1)
   
   C(1) = ANY( Mesh % Nodes % x /= x ) 
   C(2) = ANY( Mesh % Nodes % y /= y )  
   C(3) = ANY( Mesh % Nodes % z /= z )  

   ! This version is perhaps too liberal 
   Mesh % MeshDim = COUNT( C )
   Mesh % MaxDim = 0
   DO i=1,3
     IF( C(i) ) Mesh % MaxDim = i
   END DO
      
   CALL Info('SetMeshDimension','Dimension of mesh is: '//I2S(Mesh % MeshDim),Level=8)
   CALL Info('SetMeshDimension','Max dimension of mesh is: '//I2S(Mesh % MaxDim),Level=8)

 END SUBROUTINE SetMeshDimension

 
 !------------------------------------------------------------------------------
 !> Function to load mesh from disk.
 !------------------------------------------------------------------------------
 FUNCTION LoadMesh2( Model, MeshDirPar, MeshNamePar,&
     BoundariesOnly, NumProcs, MyPE, Def_Dofs, mySolver, &
     LoadOnly ) RESULT( Mesh )
   !------------------------------------------------------------------------------
   USE PElementMaps, ONLY : GetRefPElementNodes

   IMPLICIT NONE

   CHARACTER(LEN=*) :: MeshDirPar,MeshNamePar
   LOGICAL :: BoundariesOnly    
   INTEGER, OPTIONAL :: numprocs,mype,Def_Dofs(:,:), mySolver
   TYPE(Mesh_t),  POINTER :: Mesh
   TYPE(Model_t) :: Model
   LOGICAL, OPTIONAL :: LoadOnly 
   !------------------------------------------------------------------------------    
   INTEGER :: i,j,k,n
   INTEGER :: BaseNameLen, Save_Dim
   LOGICAL :: GotIt, Found
   TYPE(Element_t), POINTER :: Element
   TYPE(Matrix_t), POINTER :: Projector
   LOGICAL :: parallel, LoadNewMesh
   TYPE(ValueList_t), POINTER :: VList
   CHARACTER(:), ALLOCATABLE :: FileName
   CHARACTER(*), PARAMETER :: Caller='LoadMesh'

   Mesh => Null()
   
   n = LEN_TRIM(MeshNamePar)
   DO WHILE (MeshNamePar(n:n)==CHAR(0).OR.MeshNamePar(n:n)==' ')
     n=n-1
   END DO
   IF(NumProcs<=1) THEN
     INQUIRE( FILE=MeshNamePar(1:n)//'/mesh.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Fatal(Caller,'Requested mesh > '//MeshNamePar(1:n)//' < does not exist!')
     END IF
     CALL Info(Caller,'Loading serial mesh!',Level=8)
    
   ELSE
     INQUIRE( FILE=MeshNamePar(1:n)//'/partitioning.'// & 
         i2s(Numprocs)//'/part.1.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Warn(Caller,'Requested mesh > '//MeshNamePar(1:n)//' < in partition '&
           //I2S(MyPe)//' does not exist!')
       RETURN
     END IF
     CALL Info(Caller,'Loading parallel mesh for '//I2S(Numprocs)//' partitions',Level=8)
   END IF
     
   Parallel = .FALSE.
   IF ( PRESENT(numprocs) .AND. PRESENT(mype) ) THEN
     IF ( numprocs > 1 ) Parallel = .TRUE.
   END IF

   Mesh => AllocateMesh()

   ! Get sizes of mesh structures for allocation
   !--------------------------------------------------------------------
   CALL LoadMeshStep( 1, Mesh, MeshNamePar, mype, numprocs, Parallel )

   ! Initialize and allocate mesh structures
   !---------------------------------------------------------------------
   IF( BoundariesOnly ) Mesh % NumberOfBulkElements = 0
   CALL InitializeMesh( Mesh )

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
   FileName = MeshNamePar(1:BaseNameLen)//'/mesh.names'
   CALL ReadTargetNames( Model, FileName )

   ! Map bodies using Target Bodies and boundaries using Target Boundaries.
   ! This must be done before the element definitions are studied since
   ! then the pointer should be to the correct body index. 
   !------------------------------------------------------------------------
   CALL MapBodiesAndBCs()

   ! Read parallel mesh information: shared nodes
   !------------------------------------------------------------------
   CALL LoadMeshStep( 5 )

   ! Set default internal/external BCs. This must be after the previous load mesh
   ! since only there the shared nodes are loaded, and this info is used to decide
   ! whether a boundary element is internal or external.
   !------------------------------------------------------------------------------
   CALL MapInternalExternalBCs()


   ! If requested split quadrilaterals to triangles.
   !-----------------------------------------------------------------------
   CALL SplitMeshQuads(Mesh, Model % Simulation )
   
   ! Create new boundaries on intersection of boundaries or bodies.
   ! This way the original mesh does not need to include the BCs
   ! initially.
   !-------------------------------------------------------------------
   CALL CreateIntersectionBCs(Model, Mesh)

   ! Sometimes we need boundaries that do not exist in the original mesh.
   ! Then we may create boundaries based on some geometric rules. 
   !--------------------------------------------------------------------
   CALL TagBCsUsingRule(Model, Mesh)
   
   ! Create the discontinuous mesh that accounts for the jumps in BCs
   ! This must be created after the whole mesh has been read in and 
   ! bodies and bcs have been mapped to full operation.
   ! To consider non-nodal elements it must be done before them.
   !--------------------------------------------------------------------
   CALL CreateDiscontMesh(Model,Mesh)

   ! Deallocate some stuff no longer needed
   !------------------------------------------------------------------
   CALL LoadMeshStep( 6 )

   CALL Info(Caller,'Loading mesh done',Level=8)
   
   IF( PRESENT( LoadOnly ) ) THEN
     CALL Info(Caller,'Only loading mesh, saving final preparation for later!',Level=12)     
     IF( LoadOnly ) RETURN
   END IF

   IF( PRESENT( mySolver ) ) THEN     
     VList => Model % Solvers(mySolver) % Values
   ELSE
     VList => Model % Simulation
   END IF
   IF(.NOT. ListGetLogical( VList,'Finalize Meshes Before Extrusion',Found ) ) THEN
     ! The final preparation for the mesh (including dof definitions) will be
     ! done only after the mesh has been extruded. 
     IF( ListCheckPrefix( VList,'Extruded Mesh') ) THEN
       CALL Info(Caller,'This mesh will be extruded, skipping finalization',Level=12)
       RETURN
     END IF
   END IF
   
   CALL PrepareMesh(Model,Mesh,Parallel,Def_Dofs,mySolver)

   CALL Info(Caller,'Preparing mesh done',Level=10)

   IF( Parallel ) CALL RadiationParallelMeshDistribute(Mesh, NumProcs)
   
 CONTAINS


   !------------------------------------------------------------------------------
   ! Map bodies and boundaries as prescirbed by the 'Target Bodies' and 
   ! 'Target Boundaries' keywords.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapBodiesAndBCs()

     TYPE(Element_t), POINTER :: Element
     INTEGER, ALLOCATABLE :: IndexMap(:), TmpIndexMap(:)
     INTEGER, POINTER :: Blist(:)
     INTEGER :: id,minid,maxid,body,bndry,DefaultTargetBC, DefaultTargetBody


     ! If "target bodies" is used map the bodies accordingly
     !------------------------------------------------------
     Found = .FALSE. 
     DefaultTargetBody = 0
     DO id=1,Model % NumberOfBodies
       IF( ListCheckPresent( Model % Bodies(id) % Values,'Target Bodies') ) Found = .TRUE.
       IF(ListGetLogical( Model % Bodies(id) % Values, &
           'Default Target', GotIt)) THEN
         DefaultTargetBody = id
         Found = .TRUE.
       END IF
     END DO

     IF( DefaultTargetBody /= 0 ) THEN
       CALL Info('MapBodiesAndBCs','Default Target Body: '&
           //I2S(DefaultTargetBody),Level=8)
     END IF
     
     IF( Found ) THEN
       CALL Info('MapBodiesAndBCs','Remapping bodies',Level=8)      
       minid = HUGE( minid ) 
       maxid = -HUGE( maxid ) 
       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         id = Element % BodyId
         minid = MIN( id, minid ) 
         maxid = MAX( id, maxid )
       END DO
       IF( minid > maxid ) THEN
         CALL Fatal('MapBodiesAndBCs','Body indexes are screwed!')
       END IF
       CALL Info('MapBodiesAndBCs','Minimum initial body index: '//I2S(minid),Level=6 )
       CALL Info('MapBodiesAndBCs','Maximum initial body index: '//I2S(maxid),Level=6 )

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
               CONTINUE
             ELSE IF( IndexMap( body ) /= 0 ) THEN
               CALL Warn('MapBodiesAndBCs','Multiple bodies have same > Target Bodies < entry : '&
                   //I2S(body))
             ELSE
               IndexMap( body ) = id 
             END IF
           END DO
         ELSE
           IF( DefaultTargetBody == 0 ) THEN
             IF( IndexMap( id ) /= 0 ) THEN
               CALL Warn('MapBodiesAndBCs','Unset body already set by > Target Boundaries < : '&
                   //I2S(id) )
             ELSE 
               IndexMap( id ) = id
             END IF
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

         IF( IndexMap( id ) == 0 ) THEN
           IF( DefaultTargetBody /= 0 ) THEN
             IndexMap( id ) = DefaultTargetBody
           END IF
         END IF

         Element % BodyId = IndexMap( id ) 
       END DO

       DEALLOCATE( IndexMap )
     ELSE
       CALL Info('MapBodiesAndBCs','Skipping remapping of bodies',Level=10)      
     END IF


     IF( Mesh % NumberOfBoundaryElements == 0 ) RETURN

     ! Target boundaries are usually given so this is not conditional
     !---------------------------------------------------------------
     CALL Info('MapBodiesAndBCs','Remapping boundaries',Level=8)      
     minid = HUGE( minid ) 
     maxid = -HUGE( maxid ) 
     DO i=Mesh % NumberOfBulkElements+1,&
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       id = Element % BoundaryInfo % Constraint
       minid = MIN( id, minid ) 
       maxid = MAX( id, maxid )
     END DO


     CALL Info('MapBodiesAndBCs','Minimum initial boundary index: '//I2S(minid),Level=6 )
     CALL Info('MapBodiesAndBCs','Maximum initial boundary index: '//I2S(maxid),Level=6 )
     IF( minid > maxid ) THEN
       CALL Fatal('MapBodiesAndBCs','Boundary indexes are screwed')
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
         CALL Warn('MapBodiesAndBCs','BoundaryId exceeds range')
       ELSE IF( bndry == 0 ) THEN
         CALL Warn('MapBodiesAndBCs','BoundaryId is zero')
       ELSE
         IndexMap( bndry ) = id
       END IF
     END DO

     DefaultTargetBC = 0
     DO id=1,Model % NumberOfBCs
       IF(ListGetLogical( Model % BCs(id) % Values, &
           'Default Target', GotIt)) DefaultTargetBC = id       
       IF(ListGetLogical( Model % BCs(id) % Values, &
           'Default BC', GotIt)) DefaultTargetBC = id       
       BList => ListGetIntegerArray( Model % BCs(id) % Values, &
           'Target Boundaries', GotIt )
       IF ( Gotit ) THEN
         DO k=1,SIZE(BList)
           bndry = Blist(k)
           IF( bndry > maxid ) THEN
             CONTINUE
           ELSE IF( IndexMap( bndry ) /= 0 ) THEN
             CALL Warn('MapBodiesAndBCs','Multiple BCs have same > Target Boundaries < entry : '&
                 //I2S(bndry) )
           ELSE 
             IndexMap( bndry ) = id 
           END IF
         END DO
       ELSE
         IF (ListCheckPresent(Model % BCs(id) % Values, 'Target Nodes') .OR. &
             ListCheckPresent(Model % BCs(id) % Values, 'Target Coordinates')) &
             CYCLE
         IF (IndexMap( id ) /= 0 .AND. id == DefaultTargetBC ) THEN ! DefaultTarget has been given
           CALL Warn('MapBodiesAndBCs','Default Target is a Target Boundaries entry in > Boundary Condition < : '&
               //I2S(IndexMap(id)) )
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
       CALL Info('MapBodiesAndBCs','Default Target BC: '&
           //I2S(DefaultTargetBC),Level=8)
     END IF


     DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       n = Element % TYPE % NumberOfNodes
       bndry = Element % BoundaryInfo % Constraint 

       IF( bndry > maxid .OR. bndry < minid ) THEN
         CALL Warn('MapBodiesAndBCs','Boundary index '//I2S(bndry)&
             //' not in range: '//I2S(minid)//','//I2S(maxid) )
       END IF

       IF( IndexMap( bndry ) < 0 ) THEN
         Element % BoundaryInfo % Constraint = 0
         CYCLE

       ELSE IF( IndexMap( bndry ) == 0 ) THEN
         IF( DefaultTargetBC /= 0 ) THEN
           IndexMap( bndry ) = DefaultTargetBC
         ELSE 
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



   !------------------------------------------------------------------------------
   ! Map bodies and boundaries as prescirbed by the 'Target Bodies' and 
   ! 'Target Boundaries' keywords.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapInternalExternalBCs()

     TYPE(Element_t), POINTER :: Element
     INTEGER :: id,minid,maxid,bndry,m,&
         DefaultIntBC, DefaultExtBC, cnt, cntInt, cntExt, dim

     IF( Mesh % NumberOfBoundaryElements == 0 ) RETURN

     ! Check if default internal/external BCs given
     !------------------------------------------------------------------
     DefaultIntBC = 0
     DefaultExtBC = 0
     DO id=1,Model % NumberOfBCs
       IF(ListGetLogical( Model % BCs(id) % Values, &
           'Default Internal BC', GotIt)) DefaultIntBC = id       
       IF(ListGetLogical( Model % BCs(id) % Values, &
           'Default External BC', GotIt)) DefaultExtBC = id       
     END DO
     IF(DefaultIntBC == 0 .AND. DefaultExtBC == 0) RETURN

     IF( DefaultIntBC /= 0 ) THEN
       CALL Info('MapInternalExternalBCs','Default Internal BC: '//I2S(DefaultIntBC),Level=8)
     END IF
     IF( DefaultExtBC /= 0 ) THEN
       CALL Info('MapInternalExternalBCs','Default External BC: '//I2S(DefaultExtBC),Level=8)
     END IF


     ! And finally set internal/external BCs
     !---------------------------------------
     cntInt = 0
     cntExt = 0
     dim = Mesh % MeshDim
     DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       n = Element % TYPE % NumberOfNodes
       bndry = Element % BoundaryInfo % Constraint 
       IF(bndry /= 0) CYCLE

       ! The internal/external is defined by the number of parent.
       ! This is meaningful only for dim-1 elements. Others are ignored. 
       IF(dim == 3 ) THEN
         IF( Element % TYPE % ElementCode < 300 ) CYCLE
       ELSE IF( dim == 2 ) THEN
         IF( Element % TYPE % ElementCode < 200 ) CYCLE
       END IF
       
       cnt = 0
       IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) cnt = cnt + 1
       IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) cnt = cnt + 1

       ! In parallel we may have a invalid external BC so check that it is not
       ! really an internal one.
       IF( cnt == 1 .AND. ParEnv % PEs > 1) THEN
         m = 0
         DO j=1,n
           k = Element % NodeIndexes(j)
           IF(ASSOCIATED(Mesh % ParallelInfo % NeighbourList(k) % Neighbours ) ) THEN
             IF(SIZE(Mesh % ParallelInfo % NeighbourList(k) % Neighbours) > 1) m = m+1
           END IF
         END DO
         IF(m==n) cnt = 2
       END IF       

       IF( cnt == 2 .AND. DefaultIntBC > 0 ) THEN
         cntInt = cntInt + 1
         Element % BoundaryInfo % Constraint = DefaultIntBC
       ELSE IF( cnt == 1 .AND. DefaultExtBC > 0 ) THEN
         cntExt = cntExt + 1
         Element % BoundaryInfo % Constraint = DefaultExtBC
       END IF
     END DO

     IF( cntInt /= 0 ) THEN
       CALL Info('MapInternalExternalBCs','"Default Internal BC" count: '&
           //I2S(cntInt),Level=6)
     END IF
     IF( cntExt /= 0 ) THEN
       CALL Info('MapInternalExternalBCs','"Default External BC" count: '&
           //I2S(cntExt),Level=6)
     END IF

   END SUBROUTINE MapInternalExternalBCs
   

   !------------------------------------------------------------------------------
   ! Map and scale coordinates, and increase the size of the coordinate
   ! vectors, if requested.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapCoordinates()

     REAL(KIND=dp), POINTER CONTIG :: NodesX(:), NodesY(:), NodesZ(:)
     REAL(KIND=dp), POINTER :: Wrk(:,:)
     INTEGER, POINTER :: CoordMap(:)
     REAL(KIND=dp) :: CoordScale(3)
     INTEGER :: mesh_dim, model_dim
     
     ! Perform coordinate mapping
     !------------------------------------------------------------
     CoordMap => ListGetIntegerArray( Model % Simulation, &
         'Coordinate Mapping',GotIt )
     IF ( GotIt ) THEN
       CALL Info('MapCoordinates','Performing coordinate mapping',Level=8)

       IF ( SIZE( CoordMap ) /= 3 ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'MapCoordinates', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'MapCoordinates', Message )
       END IF

       IF ( ALL( CoordMap(1:3) /= 1 ) .OR. ALL( CoordMap(1:3) /= 2 ) .OR. ALL( CoordMap(1:3) /= 3 ) ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'MapCoordinates', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'MapCoordinates', Message )
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
     CALL SetMeshDimension( Mesh )
     
     mesh_dim = Mesh % MaxDim

     ! Scaling of coordinates
     !-----------------------------------------------------------------------------
     Wrk => ListGetConstRealArray( Model % Simulation,'Coordinate Scaling',GotIt )    
     IF( GotIt ) THEN            
       CoordScale = 1.0_dp
       DO i=1,mesh_dim
         j = MIN( i, SIZE(Wrk,1) )
         CoordScale(i) = Wrk(j,1)
       END DO
       WRITE(Message,'(A,3ES10.3)') 'Scaling coordinates:',CoordScale(1:3)
       CALL Info('MapCoordinates',Message) 
       Mesh % Nodes % x = CoordScale(1) * Mesh % Nodes % x
       IF( mesh_dim > 1 ) Mesh % Nodes % y = CoordScale(2) * Mesh % Nodes % y
       IF( mesh_dim > 2 ) Mesh % Nodes % z = CoordScale(3) * Mesh % Nodes % z
     END IF

   END SUBROUTINE MapCoordinates

 !------------------------------------------------------------------------------
 END FUNCTION LoadMesh2
 !------------------------------------------------------------------------------


 !> Prepare a clean nodal mesh as it comes after being loaded from disk.
 !> Study the non-nodal elements (face, edge, DG, and p-elements)
 !> Create parallel info for the non-nodal elements
 !> Enlarge the coordinate vectors for p-elements.
 !> Generate static projector for periodic BCS.
 !-------------------------------------------------------------------
 SUBROUTINE PrepareMesh( Model, Mesh, Parallel, Def_Dofs, mySolver )
   TYPE(Model_t) :: Model
   TYPE(Mesh_t), POINTER :: Mesh
   LOGICAL :: Parallel
   INTEGER, OPTIONAL :: Def_Dofs(:,:), mySolver
   TYPE(ValueList_t), POINTER :: Vlist

   LOGICAL :: Found, DoIt
   CHARACTER(*),PARAMETER :: Caller='PrepareMesh'      
   
   IF( PRESENT( mySolver ) ) THEN     
     VList => Model % Solvers(mySolver) % Values
   ELSE
     VList => Model % Simulation
   END IF
   
   IF( Mesh % MaxDim == 0) THEN
     CALL SetMeshDimension( Mesh )
   END IF
   Model % DIMENSION = MAX( Model % DIMENSION, Mesh % MaxDim ) 

   CALL SplitMeshQuads( Mesh, Vlist ) 
   
   IF( ListGetLogical( Vlist,'Constant Stencil', Found ) ) THEN
     CALL SetEqualElementIndeces( Mesh )
   END IF

   IF( ListGetLogical( Vlist,'Increase Element Order',Found ) ) THEN
     ! We need to follow the boundary also for the new nodes of the quadratic mesh.
     CALL FollowCurvedBoundary( Model, Mesh, .FALSE. ) 
     CALL EnlargeCoordinates( Mesh ) 
     CALL FollowCurvedBoundary( Model, Mesh, .FALSE. ) 
     CALL IncreaseElementOrder( Model, Mesh )
   END IF

   CALL NonNodalElements()

   IF( Parallel ) THEN
     CALL Info(Caller,'Generating parallel communications for the non-nodal mesh',Level=20)
     CALL ResetTimer('ParallelNonNodal')
     CALL ParallelNonNodalElements()
     CALL CheckTimer('ParallelNonNodal',Level=7,Delete=.TRUE.)
   END IF

   CALL EnlargeCoordinates( Mesh ) 

   CALL FollowCurvedBoundary( Model, Mesh, .FALSE. ) 

   
   CALL GeneratePeriodicProjectors( Model, Mesh )    
   
   IF( ListGetLogical( Vlist,'Inspect Quadratic Mesh', Found ) ) THEN
     CALL InspectQuadraticMesh( Mesh ) 
   END IF
   
   IF(ListGetLogical( Model % Simulation, 'Parallel Reduce Element Max Sizes', Found ) ) THEN
     Mesh % MaxElementDOFs  = ParallelReduction( Mesh % MaxElementDOFs,2  ) 
     Mesh % MaxElementNodes = ParallelReduction( Mesh % MaxElementNodes,2 ) 
   END IF

   DoIt = ListGetLogical( Vlist,'Inspect Mesh',Found ) .OR. &
       ListGetLogical( Vlist,'Check Mesh',Found ) 
   
   IF( InfoActive(20) .OR. DoIt .OR. ListGetLogical( Vlist,'Size Info',Found ) ) THEN
     CALL PrintMeshSize( Mesh )
   END IF

   IF(DoIt) CALL CheckMeshInfo( Mesh ) 

 CONTAINS
     

   ! Check for the non-nodal element basis
   !--------------------------------------------------------
   SUBROUTINE NonNodalElements()

     INTEGER, POINTER :: EdgeDofs(:), FaceDofs(:)
     INTEGER :: i, j, k, k2, l, s, n, DGIndex, body_id, body_id0, eq_id, solver_id, el_id, &
         mat_id
     LOGICAL :: NeedEdges, Found, FoundDef0, FoundDef, FoundEq, GotIt, MeshDeps, &
         FoundEqDefs, FoundSolverDefs(Model % NumberOfSolvers), &
         FirstOrderElements, InheritDG, Hit, Stat, &
         UpdateDefDofs(Model % NumberOfSolvers)
     TYPE(Element_t), POINTER :: Element, Parent, pParent
     TYPE(Element_t) :: DummyElement
     TYPE(ValueList_t), POINTER :: Vlist
     INTEGER :: inDOFs(10,6)
     CHARACTER(MAX_NAME_LEN) :: ElementDef0, ElementDef
     
     
     EdgeDOFs => NULL()
     CALL AllocateVector( EdgeDOFs, Mesh % NumberOfBulkElements, Caller )
     FaceDOFs => NULL()
     CALL AllocateVector( FaceDOFs, Mesh % NumberOfBulkElements, Caller )     
    
     DGIndex = 0

     IF ( PRESENT(Def_Dofs) ) THEN
       inDofs = Def_Dofs
     ELSE
       InDofs = 0
       InDofs(:,1) = 1
       InDofs(:,4) = -1       
       DO s=1,Model % NumberOfSolvers
         DO i=1,6
           DO j=1,10
             inDofs(j,i) = MAX(Indofs(j,i),MAXVAL(Model % Solvers(s) % Def_Dofs(j,:,i)))
           END DO
         END DO
       END DO
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
    MeshDeps = .FALSE.  ! The order of p-basis given with a MATC function
    FoundEqDefs = .FALSE.;  FoundSolverDefs = .FALSE.

    !
    ! As a preliminary step, check if an element definition is given 
    ! in an equation section. The more common way is to give the element
    ! definition in a solver section.
    !
    DO eq_id=1,Model % NumberOFEquations
      Vlist => Model % Equations(eq_id) % Values
      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0)
      FoundEqDefs = FoundEqDefs .OR. FoundDef0

      IF (FoundDef0) THEN
        !
        ! Check if the order of p-basis is defined by calling a special
        ! MATC function:
        !
        j = INDEX(ElementDef0,'p:')
        IF (j>0 .AND. ElementDef0(j+2:j+2)=='%') MeshDeps = .TRUE.
      ELSE
        !
        ! Check if element definitions are given for each solver separately
        ! by using a special keyword construct and tag the corresponding
        ! entries in the list of the solvers. 
        ! 
        DO Solver_id=1,Model % NumberOfSolvers
          IF (PRESENT(mySolver)) THEN
            IF ( Solver_id /= mySolver ) CYCLE
          ELSE
            ! Respect definitions given in the solver section:
            IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
          END IF

          ElementDef = ListGetString(Vlist,'Element{'//i2s(solver_id)//'}',FoundDef)
          FoundSolverDefs(Solver_id) = FoundSolverDefs(solver_id) .OR. FoundDef

          IF (FoundDef) THEN
            j = INDEX(ElementDef,'p:')
            IF (j>0 .AND. ElementDef(j+2:j+2)=='%') MeshDeps = .TRUE.
          END IF
        END DO
      END IF
    END DO

    !
    ! Tag solvers for which the element definition has been given in
    ! a solver section. The function LoadModel has already read these
    ! element definitions except for cases where the order of p-basis is
    ! defined in terms of a MATC function. The array UpdateDefDofs will
    ! show whether element definitions should be re-read.
    !
    UpdateDefDofs = .TRUE.
    DO solver_id=1,Model % NumberOfSolvers
      Vlist => Model % Solvers(solver_id) % Values

      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0)

      IF (FoundDef0) THEN
        FoundSolverDefs(Solver_id) = .TRUE.

        j = INDEX(ElementDef0,'p:')
        IF (j>0 .AND. ElementDef0(j+2:j+2)=='%') THEN
          meshdeps = .TRUE.
        ELSE
          ! Solverwise element definitions have already be read in LoadModel,
          ! indicate that re-reading is not needed here
          UpdateDefDofs(Solver_id) = .FALSE.
        END IF
      ELSE
        ! If an element definition is given in an equation section, the above code
        ! does not indicate for which solvers the definition is active, so
        ! the following array update is conditional 
        IF (.NOT. FoundEqDefs) UpdateDefDofs(Solver_id) = .FALSE.
      END IF
    END DO

    ! The basic case without the order of p-basis being defined by a MATC function:
    !
    IF (.NOT.MeshDeps) THEN
      FoundDef0 = .FALSE.
      DO body_id=1,Model % NumberOfBodies
        ElementDef0 = ' '
        Vlist => Model % Bodies(body_id) % Values
        eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
        IF ( FoundEq ) THEN
          Vlist => Model % Equations(eq_id) % Values
          IF (FoundEqDefs) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

          DO solver_id=1,Model % NumberOfSolvers

            IF(PRESENT(mySolver)) THEN
              IF ( Solver_id /= mySolver ) CYCLE
            ELSE
              IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
            END IF

            ElementDef = ListGetString(Model % Bodies(body_id) % Values, &
                'Solver '//i2s(solver_id)//': Element',FoundDef )
            IF (FoundDef) THEN
              CALL Info('NonNodalElements',&
                  'Element defined in body '//i2s(Body_Id)//' for solver '//i2s(Solver_Id), Level=7) 
              CALL Info('NonNodalElements','The element definition string is '//ElementDef, Level=7)
              CALL GetMaxDefs( Model, Mesh, DummyElement, ElementDef, solver_id, body_id, Indofs )
              CYCLE
            END IF
            
            FoundDef = .FALSE.
            IF(FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//i2s(solver_id)//'}',FoundDef)

            IF ( FoundDef ) THEN
              CALL GetMaxDefs( Model, Mesh, DummyElement, ElementDef, solver_id, body_id, Indofs )
            ELSE
              IF (UpdateDefDofs(Solver_id)) THEN
                IF (.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) &
                    ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

                CALL GetMaxDefs( Model, Mesh, DummyElement, ElementDef0, solver_id, body_id, Indofs )

                IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
              ! ELSE
              !   PRINT *, 'NO NEED TO RECREATE DEF_DOFS '
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF

     ! non-nodal elements in bulk elements
     !------------------------------------------------------------
     body_id0 = -1; FoundDef=.FALSE.; FoundEq=.FALSE.
     ElementDef = ' '
     

     ! Check whether face DOFs have been generated by "-quad_face b: ..." or
     ! "-tri_face b: ..."
     !
     NeedEdges = ANY( inDOFs(9:10,5)>0 )

     DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       body_id = Element % BodyId
       n = Element % TYPE % NumberOfNodes
       
       ! Check if the order of p-basis depends on a MATC function
       IF ( Meshdeps ) THEN
         IF ( body_id/=body_id0 ) THEN
           Vlist => Model % Bodies(body_id) % Values
           eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
           ElementDef0 = ' '
         END IF

         IF ( FoundEq ) THEN
           Vlist => Model % Equations(eq_id) % Values
           FoundDef0 = .FALSE.
           IF( FoundEqDefs.AND.body_id/=body_id0 ) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

           DO solver_id=1,Model % NumberOfSolvers
             IF(PRESENT(mySolver)) THEN
               IF ( Solver_id /= mySolver ) CYCLE
             ELSE
               IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
             END IF

             FoundDef = .FALSE.
             IF (FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//i2s(solver_id)//'}',FoundDef)

             IF ( FoundDef ) THEN
               CALL GetMaxDefs( Model, Mesh, Element, ElementDef, solver_id, body_id, Indofs )
             ELSE
               IF (UpdateDefDofs(Solver_id)) THEN
                 IF (.NOT. FoundDef0.AND.FoundSolverDefs(solver_id)) &
                     ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

                 CALL GetMaxDefs( Model, Mesh, Element, ElementDef0, solver_id, body_id, Indofs )

                 IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
               END IF
             END IF
           END DO
         END IF
         body_id0 = body_id
       END IF


       el_id = Element % TYPE % ElementCode / 100

       ! Apply the elementtypes

       Element % NDOFs = n * MAX(0,inDOFs(el_id,1)) ! The count of all nodal DOFs for the element
       EdgeDOFs(i) = MAX(0,inDOFs(el_id,2))
       FaceDOFs(i) = MAX(0,inDOFs(el_id,3))

       IF ( inDofs(el_id,4) == 0 ) THEN
         inDOFs(el_id,4) = n
       END IF

       NULLIFY( Element % DGIndexes )
       IF ( inDOFs(el_id,4) > 0 ) THEN
         CALL AllocateVector( Element % DGIndexes, inDOFs(el_id,4))
         IF( indofs(el_id,4) /= Element % TYPE % NumberOfNodes ) &
             PRINT *,'Element:',Element % TYPE % ElementCode, indofs(el_id,4)
         DO j=1,inDOFs(el_id,4)
           DGIndex = DGIndex + 1
           Element % DGIndexes(j) = DGIndex
         END DO
       END IF
       Element % DGDOFs = MAX(0,inDOFs(el_id,4))
       NeedEdges = NeedEdges .OR. ANY( inDOFs(el_id,2:4)>0 )
       
       
       ! Check if given element is a p element
       IF (FirstOrderElements .AND. inDOFs(el_id,6) > 0) THEN
         CALL AllocatePDefinitions(Element)
         NeedEdges = .TRUE.
         
         ! Calculate element bubble dofs and set element p

         Element % PDefs % P = inDOFs(el_id,6)   ! NOTE: If the order of p-basis is given by
                                                 ! a MATC function, the order is here defined
                                                 ! to be the maximum order over the element
                                                 ! processed so far. This is 
                                                 ! erroneous as the resulting p-distribution  
                                                 ! thus depends on the numbering of geometric
                                                 ! entities.
         !
         ! Try to fix the issue described in the above remark in a special case 
         ! where a single element definition is given in the equation section:
         !
         IF (FoundEqDefs .AND. Model % NumberOfSolvers > 0) THEN
           ! All solvers have the same element definition, pick one of these
           ! to set the polynomial degree:
           Element % PDefs % P = Model % Solvers(1) % Def_Dofs(el_id,Body_Id,6)
         END IF

         IF ( inDOFs(el_id,5) > 0 ) THEN
           Element % BDOFs = inDOFs(el_id,5)
         ELSE
           Element % BDOFs = getBubbleDOFs(Element, Element % PDefs % P)
         END IF

         ! All elements in actual mesh are not edges
         Element % PDefs % isEdge = .FALSE.

         ! If element is of type tetrahedron and is a p element, 
         ! do the Ainsworth & Coyle trick
         IF (Element % TYPE % ElementCode == 504) CALL ConvertToACTetra(Element)
         CALL GetRefPElementNodes( Element % Type,  Element % Type % NodeU, &
             Element % Type % NodeV, Element % Type % NodeW )
       ELSE 
         ! Clear P element definitions and set manual bubbles
         Element % PDefs => NULL()
         Element % BDOFs = MAX(0,inDOFs(el_id,5))
         ! WRITE (*,*) Element % BDOFs
       END IF

       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes,Element % TYPE % NumberOfNodes )
     END DO

     InheritDG = .FALSE.
     IF( dgindex > 0 ) THEN
       InheritDG = ListCheckPresentAnyMaterial( CurrentModel,'DG Parent Material')
     END IF
     
     ! non-nodal elements in boundary elements
     !------------------------------------------------------------
     k2 = 0
     DO i = Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('NonNodalElements','Element '//I2S(i)//' not associated!')
       END IF

       IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
         CALL Fatal('NonNodalElements','Type in Element '//I2S(i)//' not associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       el_id = Element % TYPE % ElementCode / 100
       Element % NDOFs  = n * MAX(0,inDOFs(el_id,1))
       
       !
       ! NOTE: The following depends on what dofs have been introduced
       ! by using the construct "-quad_face b: ..." and
       ! "-tri_face b: ..."
       !
       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) THEN
         IF( Element % BoundaryInfo % Left % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         j = Element % BoundaryInfo % Left % ElementIndex
         Element % BDOFs = 0

         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           IF(j<1 .OR. j>SIZE(EdgeDOFs)) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Left % BoundaryInfo)) THEN
               IF(ASSOCIATED(Element % BoundaryInfo % Left % BoundaryInfo % Left)) THEN
                 j = Element % BoundaryInfo % Left % BoundaryInfo % Left % ElementIndex
                 IF(j<1 .OR. j>SIZE(EdgeDOFs)) THEN
                   k2 = k2 + 1
                 ELSE
                   Element % BDOFs = EdgeDOFs(j)
                 END IF
               ELSE
                 k2 = k2 + 1
               END IF
             ELSE
               k2 = k2 + 1
             END IF            
           ELSE
             Element % BDOFs = EdgeDOFs(j)
           END IF
         ELSE
           IF(j<1 .OR. j>SIZE(FaceDofs)) THEN
             k2 = k2 + 1
           ELSE
             Element % BDOFs = FaceDOFs(j)
           END IF
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) THEN
         IF ( Element % BoundaryInfo % Right % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         j = Element % BoundaryInfo % Right % ElementIndex
         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           IF(j<1 .OR. j>SIZE(EdgeDOFs)) THEN
             IF(ASSOCIATED(Element % BoundaryInfo % Right % BoundaryInfo)) THEN
               IF(ASSOCIATED(Element % BoundaryInfo % Right % BoundaryInfo % Left)) THEN
                 j = Element % BoundaryInfo % Right % BoundaryInfo % Left % ElementIndex
                 IF(j<1 .OR. j>SIZE(EdgeDOFs)) THEN
                   k2 = k2 + 1
                 ELSE
                   Element % BDOFs = EdgeDOFs(j)
                 END IF
               ELSE
                 k2 = k2 + 1
               END IF
             ELSE
               k2 = k2 + 1
             END IF
           ELSE
             Element % BDOFs = EdgeDOFs(j)
           END IF
         ELSE
           IF(j<1 .OR. j>SIZE(FaceDofs)) THEN
             k2 = k2 + 1
           ELSE
             Element % BDOFs = FaceDOFs(j)
           END IF
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF

       
       ! Optionally also set DG indexes for BCs
       ! It is easy for outside boundaries, but for internal boundaries
       ! we need a flag "DG Parent Material".
       IF( InheritDG ) THEN
         IF(.NOT. ASSOCIATED( Element % DGIndexes ) ) THEN
           ALLOCATE( Element % DGIndexes(n) )
           Element % DGIndexes = 0
         END IF
         
         Hit = .TRUE.
         k = 0
         DO l=1,2        
           IF(l==1) THEN
             Parent => Element % BoundaryInfo % Left
           ELSE
             Parent => Element % BoundaryInfo % Right
           END IF
           IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
           k = k + 1
           pParent => Parent
           
           mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,&
               'Material',Found )
           IF(mat_id > 0 ) THEN           
             VList => CurrentModel % Materials(mat_id) % Values
           END IF
           IF( ASSOCIATED(Vlist) ) THEN
             Hit = ListGetLogical(Vlist,'DG Parent Material',Found )
           END IF
           IF( Hit ) EXIT
         END DO
         
         IF( k == 0 ) THEN
           CALL Fatal('NonnodalElements','Cannot define DG indexes for BC!')
         ELSE IF( k == 1 ) THEN
           Parent => pParent        
         ELSE IF(.NOT. Hit ) THEN
           CALL Fatal('NonnodalElements','Cannot define DG indexes for internal BC!')       
         END IF
         
         DO l=1,n
           DO j=1, Parent % TYPE % NumberOfNodes
             IF( Element % NodeIndexes(l) == Parent % NodeIndexes(j) ) THEN
               Element % DGIndexes(l) = Parent % DGIndexes(j)
               EXIT
             END IF
           END DO
         END DO
       END IF
       
     END DO

     IF( k2 > 0 ) THEN
       CALL Warn('NonnodalElements','Element indexes beyond face or edge table: '//I2S(k2))
     END IF
     
     
     IF ( Mesh % MaxElementDOFs <= 0 ) Mesh % MaxElementDOFs = Mesh % MaxElementNodes 

     ! Override automated "NeedEdges" if requested by the user.
     !------------------------------------------------------------------------------------
     IF(PRESENT(mySolver)) THEN
       Stat = ListGetLogical(Model % Solvers(mySolver) % Values, 'Need Edges', Found)
       IF(Found) NeedEdges = Stat
     END IF

     IF( Mesh % MeshDim == 2 ) THEN
       Stat = ListGetLogical(Model % Simulation, 'Need Edges 2D', Found)
       IF(Found) NeedEdges = Stat
     END IF

     IF( Mesh % MeshDim == 3 ) THEN
       Stat = ListGetLogical(Model % Simulation, 'Need Edges 3D', Found)
       IF(Found) NeedEdges = Stat
     END IF
     
     IF ( NeedEdges ) THEN
       CALL Info('NonNodalElements','Requested elements require creation of edges',Level=8)
       CALL SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs)
     END IF

     CALL SetMeshMaxDOFs(Mesh)

     IF( ASSOCIATED(EdgeDOFs) ) DEALLOCATE(EdgeDOFs )
     IF( ASSOCIATED(FaceDOFs) ) DEALLOCATE(FaceDOFs)

     IF( Mesh % MaxFaceDofs > 0 ) THEN
       CALL Info('NonNodalElements','Face dofs max: '//I2S(Mesh % MaxFaceDofs),Level=12)
     END IF
     IF( Mesh % MaxEdgeDofs > 0 ) THEN
       CALL Info('NonNodalElements','Edge dofs max: '//I2S(Mesh % MaxEdgeDofs),Level=12)
     END IF
     IF( Mesh % MaxElementDofs > 0 ) THEN
       CALL Info('NonNodalElements','Element dofs max: '//I2S(Mesh % MaxElementDofs),Level=12)
     END IF

     BLOCK
       LOGICAL :: DelIt
       DelIt = (Mesh % MaxFaceDofs + Mesh % MaxEdgeDofs == 0 .AND. &
           Mesh % MaxElementDOFs == Mesh % MaxElementNodes )        
       IF(DelIt) THEN
         IF(ListGetLogicalAnySolver( Model,'Discontinuous Galerkin') ) THEN
           DelIt = .FALSE.
           CALL Info('NonNodalElements','We keep the edges and faces for DG')
         END IF
       END IF
       IF(DelIt .AND. (ASSOCIATED(Mesh % Edges) .OR. ASSOCIATED(Mesh % Faces))) THEN
         CALL Info('NonNodalElements','Why the heck did we allocate the edges and faces?!')
         CALL ReleaseMeshEdgeTables( Mesh )
         CALL ReleaseMeshFaceTables( Mesh )
       END IF
     END BLOCK
     
     
   END SUBROUTINE NonNodalElements


   ! When the parallel nodal neighbours have been found 
   ! perform numbering for face and edge elements as well.
   !-------------------------------------------------------------------    
   SUBROUTINE ParallelNonNodalElements()

     INTEGER :: i,j,k,n     
     TYPE(Element_t), POINTER :: Element

     ! To be on the safe side create the parallel info if it is missing.
     IF( Mesh % NumberOfNodes > 0 ) THEN
       n = SIZE( Mesh % ParallelInfo % NeighbourList )              
       ! For unset neighbours just set the this partition to be the only owner
       DO i=1,n
         IF (.NOT.ASSOCIATED(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) THEN
           CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(i) % Neighbours,1)
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = ParEnv % mype
         END IF
       END DO
     END IF
       
     ! Create parallel numbering of faces
     CALL ResetTimer('SParFaceNumbering')
     CALL SParFaceNumbering(Mesh, .TRUE. )
     CALL CheckTimer('SParFaceNumbering',Level=7,Delete=.TRUE.)

     ! Create parallel numbering for edges
     CALL ResetTimer('SParEdgeNumbering')
     CALL SParEdgeNumbering(Mesh, .TRUE.)
     CALL CheckTimer('SParEdgeNumbering',Level=7,Delete=.TRUE.)

     ! There are mainly implemented for parallel debugging.
     ! The whole sequence is only activated when "Max Output Level >= 10". 
     IF( InfoActive(10) ) THEN
       CALL Info('ParallelNonNodalElements','Number of initial nodes: '&
           //I2S(Mesh % NumberOfNodes))
       
       CALL Info('ParallelNonNodalElements','Number of initial faces: '&
           //I2S(Mesh % NumberOfFaces))
       
       CALL Info('ParallelNonNodalElements','Number of initial edges: '&
           //I2S(Mesh % NumberOfEdges))
       
       j = 0; k = 0
       DO i=1,Mesh % NumberOfNodes
         IF( SIZE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) > 1 ) THEN
           j = j + 1
           IF( Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPe ) k = k + 1
         END IF
       END DO      
       CALL Info('ParallelNonNodalElements','Number of shared nodes: '//I2S(j))
       CALL Info('ParallelNonNodalElements','Number of owned shared nodes: '//I2S(k))
            
       IF( Mesh % NumberOfFaces > 0 ) THEN
         j = 0; k = 0 
         DO i=1,Mesh % NumberOfFaces
           IF( SIZE( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours ) > 1 ) THEN
             j = j + 1 
             IF( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours(1) == ParEnv % MyPe ) k = k + 1   
           END IF
         END DO
         CALL Info('ParallelNonNodalElements','Number of shared faces: '//I2S(j))
         CALL Info('ParallelNonNodalElements','Number of owned shared faces: '//I2S(k))

#if 0
         DO i=1,Mesh % NumberOfFaces
           IF( SIZE( Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours ) == 1 ) THEN
             BLOCK
               TYPE(Element_t), POINTER :: Face
               Face => Mesh % Faces(i)
               k = 0
               DO j=1,Face % TYPE % NumberOfNodes 
                 IF( SIZE( Mesh % ParallelInfo % NeighbourList(Face % NodeIndexes(j)) % Neighbours ) > 1 ) k = k + 1 
               END DO
               IF( k == Face % TYPE % NumberOfNodes ) THEN
                 PRINT *,'Face is shared but not listed!',ParEnv % MyPe, Mesh % NumberOfFaces,i
               END IF
             END BLOCK
           ELSE
             PRINT *,'Face is shared and listed: ',ParEnv % MyPe, Mesh % NumberOfFaces,i             
           END IF
         END DO
#endif   

       END IF
       
       IF( Mesh % NumberOfEdges > 0 ) THEN
         j = 0; k = 0
         DO i=1,Mesh % NumberOfEdges
           IF( SIZE( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours ) > 1 ) THEN
             j = j + 1
             IF( Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours(1) == ParEnv % MyPe ) k = k + 1   
           END IF
         END DO
         CALL Info('ParallelNonNodalElements','Number of shared edges: '//I2S(j))
         CALL Info('ParallelNonNodalElements','Number of owned shared edges: '//I2S(k))
       END IF
     END IF


     DO i=1,Mesh % NumberOfFaces
       Mesh % MinFaceDOFs = MIN(Mesh % MinFaceDOFs,Mesh % Faces(i) % BDOFs)
       Mesh % MaxFaceDOFs = MAX(Mesh % MaxFaceDOFs,Mesh % Faces(i) % BDOFs)
     END DO
     IF(Mesh % MinFaceDOFs > Mesh % MaxFaceDOFs) Mesh % MinFaceDOFs = Mesh % MaxFaceDOFs

     DO i=1,Mesh % NumberOfEdges
       Mesh % MinEdgeDOFs = MIN(Mesh % MinEdgeDOFs,Mesh % Edges(i) % BDOFs)
       Mesh % MaxEdgeDOFs = MAX(Mesh % MaxEdgeDOFs,Mesh % Edges(i) % BDOFs)
     END DO
     IF(Mesh % MinEdgeDOFs > Mesh % MaxEdgeDOFs) Mesh % MinEdgeDOFs = Mesh % MaxEdgeDOFs

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

   
 END SUBROUTINE PrepareMesh


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
     LOGICAL, POINTER, OPTIONAL :: LTag(:)   !< Logical tag, if used
     INTEGER, POINTER, OPTIONAL :: ITag(:)   !< Integer tag, if used
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
                 Itag(i) = Itag(k) + r_i(j)
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
   TYPE(Mesh_t), POINTER :: Mesh  
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
   LOGICAL :: Found, AlreadySet, DoIt, DoBCs, DoBodies
   INTEGER :: BodyMaps, BCMaps
   CHARACTER(*), PARAMETER :: Caller = 'ReadTargetNames'

   
   DoIt = ListGetLogical( Model % Simulation,'Use Mesh Names',Found )
   IF(DoIt) THEN   
     DoBCs = .TRUE.
     DoBodies = .TRUE.
   ELSE     
     DoBCs = .FALSE.
     DoBodies = .FALSE.
   END IF

   DoIt = ListGetLogical( Model % Simulation,'Use Mesh Body Names',Found )
   IF(Found) DoBodies = DoIt   
   DoIt = ListGetLogical( Model % Simulation,'Use Mesh Boundary Names',Found ) 
   IF(Found) DoBCs = DoIt

   IF(.NOT. (DoBodies .OR. DoBCs )) RETURN
   
   BodyMaps = 0
   BCMaps = 0

   OPEN( Unit=FileUnit, File=FileName, STATUS='old', IOSTAT=iostat )
   IF( iostat /= 0 ) THEN
     CALL Fatal(Caller,'Requested the use of entity names but this file does not exits: '//TRIM(FileName))
   END IF
   
   CALL Info(Caller,'Reading names info from file: '//TRIM(FileName),Level=10)

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
       CALL Fatal(Caller,'Could not find arguments for: '//str(i1:i2))
     END IF

     AlreadySet = .FALSE.

     DO i=1,Model % NumberOfBCs
       IF(.NOT. DoBCs) CYCLE
       Vlist => Model % BCs(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches BC '//I2S(i)
         IF( AlreadySet ) THEN
           CALL Fatal(Caller,'Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Boundaries') ) THEN
           CALL Info(Caller,'> Target Boundaries < already defined for BC '//I2S(i))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Boundaries',n,ivals(1:n))
           BodyMaps = BodyMaps + 1
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO

     DO i=1,Model % NumberOfBodies
       IF(.NOT. DoBodies) CYCLE
       Vlist => Model % Bodies(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches body '//I2S(i)
         IF( AlreadySet ) THEN
           CALL Fatal(Caller,'Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Bodies') ) THEN
           CALL Info(Caller,'> Target Bodies < already defined for Body '//I2S(i))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Bodies',n,ivals(1:n))
           BCMaps = BCMaps + 1
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO
     
     IF(.NOT. AlreadySet ) THEN
       IF( ParEnv % MyPe == 0 ) THEN
         CALL Info(Caller,'Could not map name to Body nor BC: '//name0(1:i2-i1+1), Level=20)
       END IF
     END IF

   END DO
   CLOSE(FileUnit)
      
   CALL Info(Caller,'Mapped '//I2S(BodyMaps)//' body names and '//I2S(BCMaps)//' bc names to elements!')
     
 END SUBROUTINE ReadTargetNames


!------------------------------------------------------------------------------
!> This subroutine reads elementwise input data from the file mesh.elements.data 
!> and inserts the data into the structured data variable 
!> Mesh % Elements(element_id) % PropertyData. The contents of the file should
!> be arranged as
!> 
!> element: element_id_1
!> data_set_name_1: a_1 a_2 ... a_n
!> data_set_name_2: b_1 b_2 ... b_m
!> data_set_name_3: ...
!> end
!> element: ...
!> ...
!> end
!------------------------------------------------------------------------------
  SUBROUTINE ReadElementPropertyFile(FileName,Mesh)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: FileName
     TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    CHARACTER(LEN=:), ALLOCATABLE :: str
    INTEGER :: i,j,n
    INTEGER, PARAMETER :: FileUnit = 10
    REAL(KIND=dp) :: x
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementData_t), POINTER :: PD,PD1
!------------------------------------------------------------------------------
    OPEN( Unit=FileUnit, File=FileName, STATUS='old', ERR=10 )

    ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)
    DO WHILE( ReadAndTrim(FileUnit,str) )
      READ( str(9:),*) i
      IF ( i < 0 .OR. i > Mesh % NumberOFBulkElements ) THEN
        CALL Fatal( 'ReadElementPropertyFile', 'Element id out of range.' )
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
    TYPE(Mesh_t), POINTER :: Mesh
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
    TYPE(Mesh_t), POINTER :: Mesh 
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
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, ALLOCATABLE :: CandNodes(:)
    INTEGER :: NoExt
    INTEGER, POINTER :: Inds(:)

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

      ! Compute the point in the local mesh xn \in [0,1]  
      ! and get the mesh parameter for that element from
      ! external function.
      !---------------------------------------------------
      DO i=1,n
        xn = (w(i)+w(i-1))/2.0_dp
        minhn = MIN( minhn, w(i)-w(i-1) )
        h(i) = ListGetFun( FunList, FunName, xn )
        IF( h(i) < EPSILON( h(i) ) ) THEN
          CALL Fatal('FunctionUnitDivision','Given value for h(i) was negative!')
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

      ! If the maximum error is small compared to the minimum elementsize then exit
      !-----------------------------------------------------------------------------
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
      ! Generic division given by a function
      !-----------------------------------------------------------------------
      CALL FunctionUnitDivision(w,n,FunName,ParList)
    ELSE
      ! Uniform or geometric division 
      !--------------------------------------------------------------
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

  
  FUNCTION MapExtrudedMaterial(Vlist,mat0,ilayer,EndLayer) RESULT ( mat )
    TYPE(Valuelist_t), POINTER :: Vlist
    INTEGER :: mat0, mat
    INTEGER, OPTIONAL :: ilayer
    LOGICAL, OPTIONAL :: EndLayer

    TYPE(ValueList_t), POINTER, SAVE :: PrevList
    LOGICAL, SAVE :: EndMat, ExtMat
    INTEGER, POINTER, SAVE :: ExtrudedElements(:)
    INTEGER, ALLOCATABLE, SAVE :: InvExtrudedElements(:)
    INTEGER, SAVE :: nDiv, nElems
    INTEGER :: i,j
    LOGICAL :: SetMat
    
    IF(.NOT. ASSOCIATED(PrevList,Vlist)) THEN
      IF(ALLOCATED(InvExtrudedElements))  DEALLOCATE(InvExtrudedElements)
           
      PrevList => Vlist
      EndMat = ListCheckPresent( Vlist,'Extruded Mesh End Map')
      IF(EndMat) THEN
        CALL Info('MapExtrudedMaterial','Extruded Mesh will be mapped at the ends!')
      END IF        

      ExtrudedElements => ListGetIntegerArray(Vlist,'Extruded Elements',ExtMat)                
      IF(ExtMat) THEN
        nDiv = SIZE(ExtrudedElements)
        nElems = SUM(ExtrudedElements)
        ALLOCATE(InvExtrudedElements(nElems))
        InvExtrudedElements = 0
        j = 0
        DO i=1,nDiv
          IF( ExtrudedElements(i) == 0) CYCLE
          InvExtrudedElements(j+1:j+ExtrudedElements(i)) = i
          j = j+ExtrudedElements(i)
        END DO
      ELSE
        nElems = ListGetInteger(Vlist,'Extruded Mesh Layers',UnfoundFatal=.TRUE.)
      END IF
    END IF

    mat = mat0
    IF( EndMat ) THEN
      IF(ilayer < 1 .OR. ilayer > nElems ) THEN
        CALL Fatal('MapExtrudedMaterial','Invalid body id: '//I2S(ilayer))
      END IF

      SetMat = .FALSE.
      IF( ExtMat ) THEN
        j = InvExtrudedElements(ilayer)
        SetMat = (j==1 .OR. j==nDiv)
      ELSE
        SetMat = (ilayer==1 .OR. ilayer==nElems)
      END IF             
      IF(SetMat) mat = NINT(ListGetFun( Vlist,'Extruded Mesh End Map',1.0_dp * mat0 ) )
    END IF

  END FUNCTION MapExtrudedMaterial

  
  
  SUBROUTINE CheckPointElementParents(Mesh)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Found
    INTEGER :: Misses(3)
    TYPE(Element_t), POINTER :: Element, Parent
    INTEGER :: i,j,t,t1,t2
    INTEGER, ALLOCATABLE :: OneOwner(:)
    
    t1 = Mesh % NumberOfBulkElements
    t2 = Mesh % NumberOfBoundaryElements

    Misses = 0
    DO t=t1+1,t1+t2
      Element => Mesh % Elements(t)
      IF(Element % TYPE % NumberOfNodes > 1) CYCLE
      IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) THEN
        Misses(1) = Misses(1) + 1
        CYCLE
      END IF
      Parent => Element % BoundaryInfo % Left
      IF(ASSOCIATED(Parent)) THEN
        i = Element % NodeIndexes(1)
        IF(ALL(Parent % NodeIndexes /= i)) Misses(2) = Misses(2) + 1
      END IF
      Parent => Element % BoundaryInfo % Right
      IF(ASSOCIATED(Parent)) THEN
        i = Element % NodeIndexes(1)
        IF(ALL(Parent % NodeIndexes /= i)) Misses(3) = Misses(3) + 1
      END IF      
    END DO

    i = SUM(Misses)
    IF(i == 0) RETURN

    IF( i > 0 ) THEN
      CALL Info('CheckPointElementParents',&
          'We have point '//I2S(i)//' elements with faulty parents!')
    END IF
    
    ALLOCATE(OneOwner(Mesh % NumberOfNodes))
    OneOwner = 0
    DO t=1,t1      
      Element => Mesh % Elements(t)
      DO i=1,Element % TYPE % NumberOfNodes
        j = Element % NodeIndexes(i)
        IF(OneOwner(j)==0) OneOwner(j) = t
      END DO
    END DO

    DO t=t1+1,t1+t2
      Element => Mesh % Elements(t)
      IF(Element % TYPE % NumberOfNodes > 1) CYCLE
      Parent => Element % BoundaryInfo % Left
      IF(ASSOCIATED(Parent)) THEN
        i = Element % NodeIndexes(1)
        IF(ALL(Parent % NodeIndexes /= i)) THEN
          Element % BoundaryInfo % Left => Mesh % Elements(OneOwner(i))
        END IF
      END IF
      Element % ElementIndex = t
    END DO
    
  END SUBROUTINE CheckPointElementParents

  
  ! Collect here the routines that defines the division in the exruded direction.
  !-----------------------------------------------------------------------------  
  SUBROUTINE ExtrudedDivision(Vlist, nlevels, Wtable)
    TYPE(ValueList_t), POINTER :: Vlist
    INTEGER :: nlevels
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)

    LOGICAL :: Found, GotLimits    
    REAL(KIND=dp) :: q,zmin,zmax,z
    INTEGER :: i,j,k,nDiv
    REAL(KIND=dp), POINTER :: ExtrudedLimits(:,:), ExtrudedSizes(:,:), ExtrudedRatios(:,:)
    INTEGER, POINTER :: ExtrudedElements(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtmp(:)
    
    nlevels = ListGetInteger(Vlist,'Extruded Mesh Layers',Found)
    IF( .NOT. Found ) THEN
      nlevels = ListGetInteger(Vlist,'Extruded Mesh Levels',Found)-1 
      IF(Found) THEN
        CALL ListAddNewInteger(Vlist,'Extruded Mesh Layers',nlevels)     
      END IF
    END IF
    IF(Found ) THEN
      q = ListGetCReal(Vlist,'Extruded Mesh Ratio',Found )
      IF(.NOT. Found) q = 1.0_dp
      ALLOCATE(Wtable(0:nlevels))
      CALL UnitSegmentDivision(Wtable,nlevels,Vlist)
      zmin = ListGetCReal(Vlist,'Extruded Min Coordinate',Found )
      zmax = ListGetCReal(Vlist,'Extruded Max Coordinate',Found )
      IF(.NOT. Found) zmax = 1.0_dp

      Wtable = zmin + (zmax-zmin) * Wtable      
    ELSE
      ExtrudedElements => ListGetIntegerArray(Vlist,'Extruded Elements',Found)                
      IF(.NOT. Found ) CALL Fatal('ExtrudedDivision','We should not even be here!')      
      nDiv = SIZE(ExtrudedElements) 
      
      ExtrudedLimits => ListGetConstRealArray(Vlist,'Extruded Limits',GotLimits) 
      IF(GotLimits) THEN
        IF(SIZE(ExtrudedLimits,1) /= nDiv+1 .OR. SIZE(ExtrudedLimits,2) /= 1) THEN
          CALL Fatal('ExtrudedDivision','Incompatible size for "Extruded Limits"')
        END IF
      ELSE
        ExtrudedSizes => ListGetConstRealArray(Vlist,'Extruded Sizes',Found ) 
        IF(.NOT. Found) THEN
          CALL Fatal('ExtrudedDivision','Give either "Extruded Limits" or "Extruded Sizes"!')
        END IF
        IF(SIZE(ExtrudedSizes,1) /= nDiv .OR. SIZE(ExtrudedSizes,2) /= 1) THEN
          CALL Fatal('ExtrudedDivision','Incompatible size for "Extruded Sizes"')
        END IF
      END IF

      ExtrudedRatios => ListGetConstRealArray(Vlist,'Extruded Ratios',Found)                
      IF(Found) THEN
        IF(SIZE(ExtrudedRatios,1) /= nDiv .OR. SIZE(ExtrudedRatios,2) /= 1) THEN
          CALL Fatal('ExtrudedDivision','Incompatible size for "Extruded Elements"')
        END IF
      END IF

      i = MAXVAL(ExtrudedElements)
      nlevels = SUM(ExtrudedElements)
      ALLOCATE(Wtable(0:nlevels),Wtmp(0:i))
      j = 0
      q = 1.0_dp
      DO i=1,nDiv
        IF(ASSOCIATED(ExtrudedRatios)) q = ExtrudedRatios(i,1)        

        k = ExtrudedElements(i)
        CALL GeometricUnitDivision(Wtmp,k,q)

        IF(GotLimits) THEN          
          Wtable(j:j+k) = ExtrudedLimits(i,1) + &
              Wtmp(0:k)*(ExtrudedLimits(i+1,1)-ExtrudedLimits(i,1))
        ELSE
          Wtable(j:j+k) = z + ExtrudedSizes(i,1)*Wtmp(0:k) 
          z = z + ExtrudedSizes(i,1)
        END IF
        j = j + k
      END DO
      DO i=0,nlevels
        WRITE( Message, '(A,I0,A,ES12.4)') 'w(',i,') : ',wTable(i)
        CALL Info('ExtrudedDivision', Message )
      END DO

      !CALL ListAddNewConstReal(Vlist,'Extruded Min Coordinate',Wtable(0) )
      !CALL ListAddNewConstReal(Vlist,'Extruded Max Coordinate',Wtable(nlevels) )
    END IF
    
    IF(nlevels < 2) THEN
      CALL Fatal('ExtrudedDivision','There must be at least two "Extruded Mesh Layers"!')
    END IF    
      
  END SUBROUTINE ExtrudedDivision


  ! Enable skew for extruded or initially 3D mesh, mainly intended for electrical
  ! machines. This is a library routine since we may want to perform skew right
  ! after the extrusion, if the mesh is further to be split into other elements. 
  !-----------------------------------------------------------------------------  
  SUBROUTINE SetMeshSkew(Mesh, Vlist )
    TYPE(ValueList_t), POINTER :: Vlist
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: RotorRad, AngleCoeff, RotorSkew, StatorSkew
    REAL(KIND=dp) :: zmin, zmax, Coord(3), zloc, alpha, minskew, maxskew
    LOGICAL :: Found, GotSkewFun, GotSkew, IsRotor
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    INTEGER :: NoNodes, elem, n, i, j, NodeIndex(1)
    LOGICAL :: SkewDone = .FALSE.        
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: RotorBodies(:)
    CHARACTER(*), PARAMETER :: Caller="SetMeshSkew"   

    SAVE SkewDone

    IF(SkewDone) THEN
      CALL Info(Caller,'Skew already done!',Level=10)
      RETURN
    END IF

    RotorBodies => ListGetIntegerArray( Vlist,'Rotor Bodies',Found )
    IF(.NOT. ASSOCIATED(RotorBodies) ) THEN
      RotorRad = ListGetCReal(Vlist,'Rotor Radius',Found )
      IF(.NOT. Found) THEN
        CALL Info(Caller,'Neither "Rotor Radius" or "Rotor Bodies" given!',Level=10)
        RETURN
      END IF
    END IF
    
    IF( ListGetLogical( Vlist,'Rotate in Radians',Found ) ) THEN
      CALL Info(Caller,'Using radians for skew!',Level=10)
      AngleCoeff = 1.0_dp
    ELSE
      CALL Info(Caller,'Using degrees for skew!',Level=10)
      AngleCoeff = PI / 180.0_dp
    END IF

    RotorSkew = AngleCoeff * ListGetCReal(Vlist,'Rotor Skew',GotSkew )
    GotSkewFun = ListCheckPresent( Vlist,'Rotor Skew Function')
    StatorSkew = AngleCoeff * ListGetCReal(Vlist,'Stator Skew',Found )
    GotSkew = GotSkew .OR. GotSkewFun .OR. Found
    IF(.NOT. GotSkew) THEN
      CALL Info(Caller,'No settings for skew given!',Level=10)
      RETURN
    END IF

    NoNodes = Mesh % NumberOfNodes

    zmin = ListGetCReal( Vlist,'Rotor Skew Min Coordinate',Found ) 
    IF(.NOT. Found) THEN
      zmin = ListGetCReal( Vlist,'Extruded Min Coordinate',Found ) 
    END IF
    IF(.NOT. Found) THEN
      zmin = MINVAL(Mesh % Nodes % z(1:NoNodes))
      zmin = ParallelReduction(zmin,1)
    END IF

    zmax = ListGetCReal( Vlist,'Rotor Skew Max Coordinate',Found ) 
    IF(.NOT. Found) THEN
      zmax = ListGetCReal( Vlist,'Extruded Max Coordinate',Found ) 
    END IF
    IF(.NOT. Found) THEN
      zmax = MAXVAL(Mesh % Nodes % z(1:NoNodes))
      zmax = ParallelReduction(zmax,2)
    END IF
    
    WRITE(Message,'(A,2ES12.3)') 'Coordinate range for extrusion:',zmin,zmax    
    CALL Info(Caller,Message)
    
    ALLOCATE(NodeDone(NoNodes))
    NodeDone = .FALSE.

    maxskew = -HUGE(maxskew)
    minskew = HUGE(minskew)
    
    DO elem = 1,Mesh % NumberOfBulkElements      
      Element => Mesh % Elements(elem)
      n = Element % TYPE % NumberOfNodes

      Coord(1) = SUM(Mesh % Nodes % x(Element % NodeIndexes)) / n
      Coord(2) = SUM(Mesh % Nodes % y(Element % NodeIndexes)) / n
      Coord(3) = SUM(Mesh % Nodes % z(Element % NodeIndexes)) / n

      IF(ASSOCIATED(RotorBodies)) THEN
        IsRotor = ANY( RotorBodies == Element % BodyId ) 
      ELSE
        IsRotor = (Coord(1)**2+Coord(2)**2 < RotorRad**2) 
      END IF

      DO i=1,n
        j = Element % NodeIndexes(i)
        NodeIndex(1) = j
        IF(.NOT. NodeDone(j)) THEN
          Coord(1) = Mesh % Nodes % x(j)
          Coord(2) = Mesh % Nodes % y(j)
          Coord(3) = Mesh % Nodes % z(j)

          ! Skew is not constant, perform it for each node 1st if requested. 
          zloc = (coord(3)-zmin)/(zmax-zmin)

          ! By construction this must be in [0,1]
          zloc = MAX(0.0_dp,MIN(1.0_dp,zloc))

          IF( IsRotor ) THEN
            IF(GotSkewFun) THEN
              alpha = AngleCoeff * ListGetFun( Vlist,'Rotor Skew Function',zloc)                
            ELSE
              alpha = (zloc-0.5_dp) * RotorSkew
            END IF
          ELSE
            alpha = (zloc-0.5_dp) * StatorSkew 
          END IF

          maxskew = MAX(alpha, maxskew)
          minskew = MIN(alpha, minskew)
          
          Mesh % Nodes % x(j) = Coord(1)*COS(alpha) - Coord(2)*SIN(alpha)
          Mesh % Nodes % y(j) = Coord(1)*SIN(alpha) + Coord(2)*COS(alpha)        
          NodeDone(j) = .TRUE.
        END IF
      END DO
    END DO

    SkewDone = .TRUE.

    IF(InfoActive(10)) THEN
      IF(GotSkewFun) THEN
        minskew = (180.0/PI) * ParallelReduction(minskew,1)
        maxskew = (180.0/PI) * ParallelReduction(maxskew,2)
        WRITE(Message,'(A,2ES12.3)') 'Rotor skew done with range (degrees): ',minskew,maxskew
      ELSE
        WRITE(Message,'(A,2ES12.3)') 'Rotor skew done with total angle: ',(180.0/PI)*RotorSkew
      END IF
      CALL Info(Caller,Message)
    END IF

  END SUBROUTINE SetMeshSkew
    
  
!------------------------------------------------------------------------------
!> Given a 2D mesh extrude it to be 3D. The 3rd coordinate will always
!> be at the interval [0,1]. Therefore the adaptation for different shapes
!> must be done with StructuredMeshMapper, or some similar utility. 
!> The top and bottom surface will be assigned Boundary Condition tags
!> with indexes one larger than the maximum used on by the 2D mesh. 
!> NOTE: This function handles NDOFs of the element structure in a way
!>       which is not consistent with "Element = n:N ...", with N>1 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrude(Mesh_in, Vlist) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh_in, Mesh_out
    TYPE(ValueList_t), POINTER :: Vlist
!------------------------------------------------------------------------------
    CHARACTER(:), ALLOCATABLE :: ExtrudedMeshName
    INTEGER :: i,j,k,l,n,cnt,ind(8),max_baseline_bid,max_bid,l_n,max_body,&
        ExtrudedCoord,dg_n,totalnumberofelements
    TYPE(Element_t), POINTER :: Elem_in, Elem_out
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    INTEGER :: nnodes,gnodes,gelements,ierr,bcignored,cnt101
    LOGICAL :: isParallel, Found, PreserveBaseline, Rotational, Rotate2Pi, CollectExtrudedBCs
    REAL(KIND=dp)::CurrCoord 
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
    INTEGER, POINTER :: BCLayers(:), TmpLayers(:)
    INTEGER :: NoBCLayers, bcoffset, baseline0, bclevel, BaseLineLayer, bcind, &
        m, max_bid0, in_levels, nlev
    INTEGER :: BcCounter(100)    
    LOGICAL :: GotBCLayers, DoCount
    CHARACTER(*), PARAMETER :: Caller="MeshExtrude"   

!------------------------------------------------------------------------------

    Mesh_out => AllocateMesh()
    isParallel = ( ParEnv % PEs > 1 )

    ! Create the division for the 1D mesh
    !--------------------------------------------
    CALL ExtrudedDivision(Vlist,nlev,Wtable)    
    CALL Info(Caller,'Extruding '//I2S(nlev)//' element layers on: '//TRIM(Mesh_in % Name),Level=10)
    in_levels = nlev-1
        
    ! Generate volume nodal points:
    ! -----------------------------
    n = Mesh_in % NumberOfNodes
    nnodes = (in_levels+2)*n
    gnodes = nnodes

    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )

    gelements = Mesh_in % NumberOfBulkElements

    ! There are some meshes with corrupted owners for 101 elements!
    ! This checks these nodes. 
    CALL CheckPointElementParents(Mesh_in)
    
    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo
      
      IF(.NOT. ASSOCIATED( PI_in ) ) CALL Fatal(Caller,'PI_in not associated!')
      IF(.NOT. ASSOCIATED( PI_out ) ) CALL Fatal(Caller,'PI_out not associated!')
            
      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % GInterface(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      IF(.NOT. ASSOCIATED( PI_in % NeighbourList ) ) THEN
        CALL Fatal(Caller,'Neighnours not associated!')
      END IF

      ! For unset neighbours just set the this partition to be the only owner
      DO i=1,Mesh_in % NumberOfNodes
        IF (.NOT.ASSOCIATED(PI_in % NeighbourList(i) % Neighbours)) THEN
          CALL AllocateVector(PI_in % NeighbourList(i) % Neighbours,1)
          PI_in % NeighbourList(i) % Neighbours(1) = ParEnv % Mype
        END IF
      END DO
          
      j=0
      DO i=1,Mesh_in % NumberOfNodes
        IF (PI_in % NeighbourList(i) % &
            Neighbours(1) == ParEnv % MyPE ) j=j+1
      END DO

      CALL MPI_ALLREDUCE(j,gnodes,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
      
      j=0
      DO i=1,Mesh_in % NumberOfBulkElements
        IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
      END DO
      
      CALL MPI_ALLREDUCE(j,gelements,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF

    CALL Info(Caller,'Global count of original elements: '//I2S(gelements),Level=12)
    CALL Info(Caller,'Number of nodes for extruded mesh: '//I2S(nnodes),Level=12)

    ExtrudedCoord = ListGetInteger( CurrentModel % Simulation,'Extruded Coordinate Index', &
        Found, minv=1,maxv=3 )
    IF(.NOT. Found) ExtrudedCoord = MIN(3,Mesh_in % MeshDim + 1)
    CALL Info(Caller,'Extrusion in direction of dimension: '//I2S(ExtrudedCoord),Level=12)
    
    IF( ExtrudedCoord == 1 ) THEN
      ActiveCoord => Mesh_out % Nodes % x
    ELSE IF( ExtrudedCoord == 2 ) THEN
      ActiveCoord => Mesh_out % Nodes % y
    ELSE IF( ExtrudedCoord == 3 ) THEN
      ActiveCoord => Mesh_out % Nodes % z
    END IF

    PreserveBaseline = ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found )

    CollectExtrudedBCs = ListGetLogical( CurrentModel % Simulation,'Extruded BCs Collect',Found )
    
    Rotate2Pi = .FALSE.
    Rotational = ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Rotational',Found )    
    IF( Rotational ) THEN
      Rotate2Pi = ( ABS(MAXVAL(Wtable)-MINVAL(Wtable) - 2*PI) < 1.0d-3*PI )
      IF( Rotate2Pi ) CALL Info(Caller,'Perfoming full 2Pi rotation',Level=6)
    END IF

    ! This sets the BC layers.
    ! We honor the old way of assuming just bottom and top layer so the internal BCs are
    ! set as additional layers between.
    TmpLayers => ListGetIntegerArray( CurrentModel % Simulation,'Extruded BC Layers', GotBCLayers ) 
    IF( GotBCLayers ) THEN
      NoBCLayers = 2 + SIZE( TmpLayers )      
    ELSE
      NoBCLayers = 2
    END IF
    ALLOCATE(BCLayers(NoBCLayers))
    BCLayers(1) = 0
    BCLayers(NoBCLayers) = in_levels+1    
    IF( GotBCLayers ) THEN
      CALL Info(Caller,'There will be total of '//I2S(NoBCLayers)//' layers with BCs',Level=8)
      BCLayers(2:NoBCLayers-1) = TmpLayers
      DO i=1,NoBCLayers-1
        IF(BCLayers(i) >= BCLayers(i+1)) THEN
          CALL Fatal(Caller,'BC layers should be in increasing order')
        END IF
      END DO
    END IF    

    DoCount = .FALSE.
    BCCounter = 0 
    BaseLineLayer = 0
    baseline0 = 0
    IF( PreserveBaseline ) THEN
      BaseLineLayer = ListGetInteger( CurrentModel % Simulation,'Extruded Baseline Layer', Found )
      IF(.NOT. Found) BaseLineLayer = 1
      IF( BaseLineLayer > NoBCLayers ) THEN
        CALL Fatal(Caller,"'Extruded Baseline Layer' cannot exceed: "//I2S(NoBCLayers)) 
      END IF
      CALL Info(Caller,'Baseline will be set to layer '//I2S(BaselineLayer),Level=8)
    END IF
    
    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Maximum body index in original mesh: '//I2S(max_body),Level=6)

    max_bid0 = 0
    DO j=1,Mesh_in % NumberOfBoundaryElements
      k = j + Mesh_in % NumberOfBulkElements
      Elem_in => Mesh_in % Elements(k)
      IF(.NOT. ASSOCIATED(Elem_in % BoundaryInfo)) CYCLE
      bcind = Elem_in % BoundaryInfo % constraint 
      max_bid0 = MAX(max_bid0,bcind)
    END DO        
    IF(isParallel) THEN
      j = max_bid0
      CALL MPI_ALLREDUCE(j,max_bid0,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Maximum boundary index in original mesh: '//I2S(max_bid0),Level=6)
    

    ! Create the nodes (and in parallel their global indexes).
    ! This assumes exacyly same distribution for each extruded node. 
    cnt=0
    DO i=0,in_levels+1

      ! If we rotate full 2Pi then we have natural closure!
      IF( Rotate2Pi ) THEN
        IF( i == in_levels+1) EXIT
      END IF      
      CurrCoord = Wtable( i ) 
      
      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          PI_out % GInterface(cnt) = PI_in % GInterface(j)

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(&
               SIZE(PI_in % NeighbourList(j) % Neighbours)))
          PI_out % NeighbourList(cnt) % Neighbours = &
            PI_in % NeighbourList(j) % Neighbours
          PI_out % GlobalDOFs(cnt) = PI_in % GlobalDOFs(j)+i*gnodes
        END IF

      END DO
    END DO
    Mesh_out % NumberOfNodes = cnt
    Mesh_out % Nodes % NumberOfNodes = cnt

    ! For rotational geometry map the coordinates. 
    IF( Rotational ) THEN
      BLOCK
        REAL(KIND=DP) :: x,y,z,r        
        DO i=1,cnt          
          x = Mesh_out % Nodes % x(i)
          y = Mesh_out % Nodes % y(i)
          z = Mesh_out % Nodes % z(i)

          Mesh_out % Nodes % x(i) = COS(z) * x
          Mesh_out % Nodes % y(i) = SIN(z) * x
          Mesh_out % Nodes % z(i) = y
        END DO
      END BLOCK
    END IF
        
    ! Warn about 101 elements:
    ! -------------------------
    cnt101 = 0
    bcignored = 0
    DO i=Mesh_in % NumberOfBulkElements+1, &
        Mesh_in % NumberOfBulkElements+Mesh_in % NumberOfBoundaryElements
      Elem_in => Mesh_in % Elements(i)
      IF(Elem_in % TYPE % ElementCode == 101) cnt101 = cnt101 + 1
      IF(Elem_in % BoundaryInfo % Constraint == 0 ) bcignored = bcignored + 1
    END DO

    IF(isParallel) THEN
      j=cnt101
      CALL MPI_ALLREDUCE(j,cnt101,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
      j=bcignored
      CALL MPI_ALLREDUCE(j,bcignored,1,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF
       
    IF( bcignored > 0 ) THEN
      CALL Info(Caller,"WARNING: We are skipping '//I2S(bcignored)//&
          ' non-defined BC elements in extrusion!",Level=3)
    END IF       
    IF( cnt101 > 0 ) THEN
      CALL Info(Caller,"WARNING: Historically 101's were extruded as is, now they become 202's!",Level=3)
    END IF
    
    ! Compute total number of elements needed
    ! extruded bulk + extruded bc elements
    n = Mesh_in % NumberOfBulkElements + Mesh_in % NumberOfBoundaryElements
    totalnumberofelements = n*(in_levels+1) 
    IF(.NOT. Rotate2Pi ) THEN
      ! new layer bc's
      totalnumberofelements = totalnumberofelements + NoBCLayers * Mesh_in % NumberOfBulkElements
    END IF
    IF (PreserveBaseline) THEN
      ! additional baseline elements, if requested
      totalnumberofelements = totalnumberofelements + Mesh_in % NumberOfBoundaryElements
    END IF
    ALLOCATE(Mesh_out % Elements(totalnumberofelements))
    
    ! Initialize all elements to zero
    DO i = 1, totalnumberofelements
      Elem_out => Mesh_out % Elements(i)      
      Elem_out % DGDOFs = 0
      Elem_out % NDOFs = 0
      Elem_out % BodyId = 0
      Elem_out % DGIndexes => NULL()
      Elem_out % PDefs => NULL()
      Elem_out % EdgeIndexes => NULL()
      Elem_out % FaceIndexes => NULL()
      Elem_out % BubbleIndexes => NULL()
    END DO
    Mesh_out % MaxElementNodes = 0

    
    ! Generate volume bulk elements:
    ! ------------------------------
    n = Mesh_in % NumberOfNodes
    cnt=0
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBulkElements
        cnt = cnt+1
        Elem_in => Mesh_in % Elements(j)
        Elem_out => Mesh_out % Elements(cnt)

        !Elem_out % BodyId = Elem_in % BodyId
        Elem_out % BodyId = MapExtrudedMaterial(Vlist,Elem_in % BodyId,i+1)        

        Elem_out % PartIndex = Elem_in % PartIndex
               
        ! If we have internal BC layers then find the correct index for the body
        IF( NoBCLayers > 2 ) THEN
          DO k=1,NoBCLayers-1
            IF(i < BCLayers(k+1) ) EXIT
          END DO
          Elem_out % BodyId = Elem_out % BodyId + max_body*(k-1)
        END IF
        
        m = Elem_in % TYPE % NumberOfNodes
        ind(1:m) = Elem_in % NodeIndexes(1:m) + i*n

        IF( Rotate2Pi .AND. i==in_levels ) THEN
          ind(m+1:2*m) = Elem_in % NodeIndexes(1:m)
        ELSE
          ind(m+1:2*m) = Elem_in % NodeIndexes(1:m)+(i+1)*n
        END IF
        m = 2*m
                
        Elem_out % NDOFs = m
        Mesh_out % MaxElementNodes = MAX(Mesh_out % MaxElementNodes,m)

        SELECT CASE(m)
        CASE(4)
          Elem_out % TYPE => GetElementType(404)
          ! We need to reorder for the quad element!
          k = ind(3); ind(3)=ind(4); ind(4) = k
        CASE(6)
          Elem_out % TYPE => GetElementType(706)
        CASE(8)
          Elem_out % TYPE => GetElementType(808)
        END SELECT

        Elem_out % GElementIndex = Elem_in % GelementIndex + gelements*i

        Elem_out % ElementIndex = cnt
        ALLOCATE(Elem_out % NodeIndexes(m)) 
        Elem_out % NodeIndexes = ind(1:m)
      END DO
    END DO
    Mesh_out % NumberOfBulkElements = cnt
    CALL Info(Caller,'Number of extruded bulk elements: '//I2S(cnt),Level=8)
      
    ! Add side boundaries with the bottom mesh boundary id's:
    ! (or shift ids if preserving the baseline boundary)
    ! -------------------------------------------------------
    max_bid = 0
    bcoffset = 0
    IF( PreserveBaseline ) THEN
      CALL Info(Caller,'Preserving original '//I2S(max_bid0)//' BCs',Level=8)
      bcoffset = max_bid0
    END IF

    CALL Info(Caller,'First extruded boundary element index: '//I2S(cnt+1),Level=20)
    
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        Elem_in => Mesh_in % Elements(k)
        bcind = Elem_in % BoundaryInfo % constraint

        ! Do not include BCs that are originally not activated
        IF(bcind==0) CYCLE

        cnt = cnt+1
        Elem_out => Mesh_out % Elements(cnt)  
        
        Elem_out = Elem_in

        Elem_out % ElementIndex = cnt        
        ALLOCATE(Elem_out % BoundaryInfo)
        Elem_out % BoundaryInfo = Elem_in % BoundaryInfo
        Elem_out % PartIndex = Elem_in % PartIndex
        
        ! Offset from possible baseline         
        bcind = bcind + bcoffset

        ! If we have internal BC layers then find the correct index for the body
        IF( NoBCLayers > 2 ) THEN
          DO k=1,NoBCLayers-1
            IF(i < BCLayers(k+1) ) EXIT
          END DO
          bcind = bcind + max_bid0*(k-1)
        END IF
        IF(DoCount .AND. bcind <= 100) BcCounter(bcind) = BcCounter(bcind) + 1
        
        Elem_out % BoundaryInfo % constraint = bcind
        max_bid = MAX(max_bid,bcind )

        m = Elem_in % TYPE % ElementCode / 100
        IF(m == 2) THEN
          Elem_out % NDOFs = 4
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(4)) 

          ind(1) = Elem_in % NodeIndexes(1)+i*n
          ind(2) = Elem_in % NodeIndexes(2)+i*n
          IF( Rotate2Pi .AND. i==in_levels ) THEN
            ind(3) = Elem_in % NodeIndexes(2)
            ind(4) = Elem_in % NodeIndexes(1)
          ELSE
            ind(3) = Elem_in % NodeIndexes(2)+(i+1)*n
            ind(4) = Elem_in % NodeIndexes(1)+(i+1)*n
          END IF
          
          Elem_out % NodeIndexes = ind(1:4)
          Elem_out % TYPE => GetElementType(404)
        ELSE IF(m == 1) THEN
          Elem_out % NDOFs = 2
          ALLOCATE(Elem_out % NodeIndexes(2))
          
          ind(1) = Elem_in % NodeIndexes(1)+i*n
          ind(2) = Elem_in % NodeIndexes(1)+(i+1)*n

          Elem_out % NodeIndexes = ind(1:2)
          Elem_out % TYPE => GetElementType(202)
        ELSE
          CALL Fatal(Caller,'Invalid number of nodes: '//I2S(m))
        END IF

        IF( bcind <= CurrentModel % NumberOfBCs) THEN
          k = ListGetInteger(CurrentModel % BCs(bcind) % Values,'Body Id',Found)
          IF(Found) Elem_out % BodyId = k
        END IF

        IF(ASSOCIATED(Elem_in % BoundaryInfo % Left)) THEN
          l = Elem_in % BoundaryInfo % Left % ElementIndex
          Elem_out % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Elem_in % BoundaryInfo % Right)) THEN
          l = Elem_in % BoundaryInfo % Right % ElementIndex
          Elem_out % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        ! Just check that we have correct parents. We had some issues here with
        ! corrupted initial meshes.
        BLOCK
          INTEGER :: ii,jj
          IF(ASSOCIATED(Elem_in % BoundaryInfo % Left)) THEN
            DO ii = 1, Elem_out % TYPE % NumberOfNodes
              jj = Elem_out % NodeIndexes(ii)
              IF( ALL( Elem_out % BoundaryInfo % Left % NodeIndexes /= jj ) ) THEN
                CALL Warn(Caller,'Node not available in left parent!')
              END IF
            END DO
          END IF
          IF(ASSOCIATED(Elem_in % BoundaryInfo % Right)) THEN
            DO ii = 1, Elem_out % TYPE % NumberOfNodes
              jj = Elem_out % NodeIndexes(ii)
              IF( ALL( Elem_out % BoundaryInfo % Right % NodeIndexes /= jj ) ) THEN
                CALL Warn(Caller,'Node not available in right parent!')
              END IF
            END DO
          END IF
        END BLOCK

        
      END DO
    END DO
        
    IF(isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Largest bc index after extruded BCs: '//I2S(max_bid),Level=8)
    IF(DoCount) PRINT *,'BCInd1:',BcCounter(1:20)
    

    ! Add start and finish planes except if we have a full rotational symmetry
    IF(Rotate2Pi ) GOTO 100 
    
    ! Add bottom, top, and possible mid boundaries:
    ! ---------------------------------------------
    CALL Info(Caller,'First plane boundary element index: '//I2S(cnt+1),Level=20)
    bcoffset = max_bid
    DO k=1,NoBCLayers

      bclevel = BCLayers(k)

      IF( PreserveBaseline ) THEN
        ! Register the starting point for parents of baseline elements
        IF(k == BaselineLayer ) THEN
          baseline0 = cnt
          CALL Info(Caller,'Baseline elements parents start from element index: '//I2S(cnt),Level=8)
        END IF
      END IF

      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        
        Elem_in => Mesh_in % Elements(i)
        Elem_out => Mesh_out % Elements(cnt)

        Elem_out = Elem_in
        Elem_out % PartIndex = Elem_in % PartIndex
        
        m = Elem_in % TYPE % NumberOfNodes
        Elem_out % NDOFs = m
        ALLOCATE(Elem_out % NodeIndexes(m))        
        ALLOCATE(Elem_out % BoundaryInfo)
        
        Elem_out % BoundaryInfo % Right => NULL()
        IF( bclevel == in_levels+1 ) THEN
          Elem_out % BoundaryInfo % Left => &
              Mesh_out % Elements((bclevel-1) * Mesh_in % NumberOfBulkElements+i)          
        ELSE
          Elem_out % BoundaryInfo % Left => &
              Mesh_out % Elements(bclevel * Mesh_in % NumberOfBulkElements+i)
          IF(bclevel > 0 ) THEN
            ! for internal BCs add the 2nd parent also!
            Elem_out % BoundaryInfo % Right => &
                Mesh_out % Elements((bclevel-1) * Mesh_in % NumberOfBulkElements+i)
          END IF
        END IF

        IF( CollectExtrudedBCs ) THEN
          bcind = bcoffset + k
        ELSE
          bcind = bcoffset + (k-1)*max_body + Elem_in % BodyId
        END IF
        IF(DoCount .AND. bcind <= 100) BcCounter(bcind) = BcCounter(bcind) + 1
        
        max_bid = MAX(max_bid,bcind )
        
        Elem_out % BoundaryInfo % Constraint = bcind        
        Elem_out % BodyId = 0

        IF( bcind <= CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcind) % Values,'Body Id',Found)
          IF(Found) Elem_out % BodyId = j
        END IF

        Elem_out % NodeIndexes = Elem_in % NodeIndexes + bclevel * n  

        Elem_out % ElementIndex = cnt
        Elem_out % TYPE => Elem_in % TYPE
      END DO
    END DO

    IF(isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1,MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Largest bc index after layer BCs: '//I2S(max_bid),Level=8)
    IF(DoCount) PRINT *,'BCInd2:',BcCounter(1:20)

    
    ! If baseline preservation is requested, these will be
    ! available in the given layer with original bc tags.
    ! We do this at the end but still use the smallest (original)
    ! bc constraint indeces here.
    ! -------------------------------------------------------
    CALL Info(Caller,'First plane boundary element index: '//I2S(cnt+1),Level=20)
    IF (PreserveBaseline ) THEN
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements
        
        Elem_In => Mesh_in % Elements(k)
        bcind = Elem_In % BoundaryInfo % Constraint 
        IF(bcind==0) CYCLE
        IF(DoCount .AND. bcind <= 100) BcCounter(bcind) = BcCounter(bcind) + 1
        
        cnt = cnt+1
        Elem_out => Mesh_out % Elements(cnt) 
        
        ALLOCATE(Elem_out % BoundaryInfo)
        Elem_out % BoundaryInfo = Elem_In % BoundaryInfo        
        Elem_out % BoundaryInfo % Constraint = bcind
        Elem_out % PartIndex = Elem_in % PartIndex
        
        Elem_out % TYPE => Elem_In % TYPE
        m = Elem_out % TYPE % ElementCode / 100
        Elem_out % NDOFs = m 

        k = BCLayers(BaselineLayer) * Mesh_in % NumberOfNodes
        ind(1:m) = Elem_in % NodeIndexes(1:m) + k
        
        ALLOCATE(Elem_out % NodeIndexes(m)) 
        Elem_out % NodeIndexes(1:m) = ind(1:m)

        Elem_out % ElementIndex = cnt
        
        IF(ASSOCIATED(Elem_In % BoundaryInfo % Left)) THEN
          l = Elem_in % BoundaryInfo % Left % ElementIndex + baseline0
          Elem_out % BoundaryInfo % Left => Mesh_out % Elements(l)
        END IF
        IF(ASSOCIATED(Elem_In % BoundaryInfo % Right)) THEN
          l = Elem_in % BoundaryInfo % Right % ElementIndex + baseline0
          Elem_out % BoundaryInfo % Right => Mesh_out % Elements(l)
        END IF
      END DO

      CALL Info(Caller,'Original baseline given by BCs: '//I2S(max_bid0))
    END IF
    IF(DoCount) PRINT *,'BCInd3:',BcCounter(1:20)
    CALL Info(Caller,'Last boundary element index: '//I2S(cnt),Level=20)

    DO i=1,cnt
      Elem_out => Mesh_out % Elements(i)
      IF( Elem_out % ElementIndex /= i) PRINT *,'mismatch: ',i,Elem_out % ElementIndex
    END DO

    
    
100 Mesh_out % NumberOfBoundaryElements = cnt-Mesh_out % NumberOfBulkElements
    
    Mesh_out % Name = Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs  = Mesh_out % MaxElementNodes
    Mesh_out % Stabilize = Mesh_in % Stabilize

    Mesh_out % MeshDim = MIN(3, Mesh_in % MeshDim + 1)
    CurrentModel % Dimension = MIN( CurrentModel % Dimension+1, 3 )
   
    DEALLOCATE( BCLayers ) 

    ! Check whether the *.sif file has included enough BCs.
    ! If not then add some for convenience.
    j = 0
    DO i=Mesh_out % NumberOfBulkElements+1, &
        Mesh_out % NumberOfBulkElements+Mesh_out % NumberOfBoundaryElements
      Elem_out => Mesh_out % Elements(i)      
      bcind = Elem_out % BoundaryInfo % Constraint
      IF(bcind==0) CYCLE
      j = MAX(bcind,j)
    END DO
    CALL Info(Caller,'Maximum bc constraint in extruded mesh: '//I2S(j))
    IF( j > CurrentModel % NumberOfBCs ) THEN
      CALL AppendMissingBCs(CurrentModel,j)
    END IF

    CALL SetMeshSkew(Mesh_out, CurrentModel % Simulation )
    
    CALL PrepareMesh( CurrentModel, Mesh_out, isParallel )
    
    ExtrudedMeshName = ListGetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
    IF(Found) THEN
      IF( ParEnv % PEs > 1 ) THEN
        ! Or WriteMeshToDiskPartitioned ? 
        CALL WriteMeshToDisk2( CurrentModel, Mesh_out, ExtrudedMeshName, ParEnv % MyPe )
      ELSE        
        CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
      END IF
    END IF
    
  CONTAINS

    SUBROUTINE AppendMissingBCs(Model,maxbc)
       TYPE(Model_t) :: Model
       INTEGER :: maxbc
       
       INTEGER :: i, NoBCs, tag 
       TYPE(BoundaryConditionArray_t), POINTER :: OldBCs(:) => NULL()
       
       NoBcs = Model % NumberOfBCs
       IF(NoBCs >= maxbc ) RETURN

       CALL Info(Caller,'Generating '//I2S(maxbc-NoBCs)//' dummy list BCs for convenience!',Level=5)
       
       OldBCs => Model % BCs(:)

       NULLIFY( Model % BCs )
       ALLOCATE( Model % BCs(maxbc) )
       
       DO i=1,NoBCs
         Model % BCs(i) % Values => OldBCs(i) % Values        
         tag = OldBCs(i) % Tag
         IF(tag == 0) tag = i 
         Model % BCs(i) % Tag = tag
       END DO
       IF (ASSOCIATED(OldBCs) .AND. NoBCs > 0) DEALLOCATE( OldBCs ) 
       DO i=NoBCs+1,maxbc
         Model % BCs(i) % Tag = i
       END DO
       DO i=1,maxbc
         IF(.NOT.ASSOCIATED(Model % BCs(i) % Values) .OR. i > NoBCs) THEN
           Model % BCs(i) % Values => ListAllocate()
           CALL ListAddString( Model % BCs(i) % Values,'Name','BC'//I2S(i))
         END IF
       END DO
       Model % NumberOfBCs = maxbc
       
     END SUBROUTINE AppendMissingBCs
    
!------------------------------------------------------------------------------
  END FUNCTION MeshExtrude
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> As the previous one except the extrusion is done in parallel for single meshes
!> that each take an internal in the extruded direction. This affects the coordinates
!> but also the communication pattern. A separate routine was made in order to avoid
!> introducing of bugs as the internal extrusion is a widely used feature. 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrudeSlices(Mesh_in, Vlist) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh_in, Mesh_out
    TYPE(ValueList_t), POINTER :: Vlist
!------------------------------------------------------------------------------
    CHARACTER(:), ALLOCATABLE :: ExtrudedMeshName
    INTEGER :: i,j,k,l,n,m,cnt,ind(8),bid,max_bid,l_n,max_body,bcid,&
        ExtrudedCoord,dg_n,totalnumberofelements,lastbc
    INTEGER, POINTER :: pInds(:)
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    TYPE(Element_t), POINTER :: Element
    INTEGER :: nnodes,gnodes,gelements,ierr,nlev,ilev,&
        nParMesh,nParExt,OrigPart,ElemCode,bodyid, newbcs
    LOGICAL :: isParallel, SingleIn, Found, TopBC, BotBC, &
        CollectExtrudedBCs, SeparateSlices, CreateInternalBCs, InternalBC
    INTEGER,ALLOCATABLE :: ChildBCs(:)
    REAL(KIND=dp)::CurrCoord 
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
    CHARACTER(*), PARAMETER :: Caller="MeshExtrudeSlices"   
!------------------------------------------------------------------------------

    ! The historical choice in_levels in annoying when we want to split the divisions.
    
    IF( ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found ) ) &
        CALL Fatal(Caller,'The slice version cannot handle "Preserve Baseline"!')
    
    IF( ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Rotational',Found ) ) &
        CALL Fatal(Caller,'The slice version cannot handle "Extruded Mesh Rotational"!')    
    
    isParallel = ( ParEnv % PEs > 1 )
    SingleIn = Mesh_in % SingleMesh
    
    ! Create the division for the 1D mesh
    !--------------------------------------------
    CALL ExtrudedDivision(Vlist,nlev,Wtable)    
    CALL Info(Caller,'Extruding '//I2S(nlev)//' element layers on: '//TRIM(Mesh_in % Name),Level=10)

    
    ! In parallel let us pick only our own share of the
    ! division. This logic makes it possible to have nonuniform divisions easily.
    ! The number of element layers is evenly distributed among partitions. 
    !-------------------------------------------------------------------------------
    IF( isParallel ) THEN
      nParExt = ParEnv % PEs 
      nParMesh = ListGetInteger( CurrentModel % Simulation,'Parallel Mesh Modulo',Found)
      IF(.NOT. Found) THEN
        nParMesh = 1
        IF(.NOT. SingleIn ) THEN
          CALL Fatal(Caller,'This routine expects either Mesh Modulo or Single Mesh!')
        END IF
      END IF
      
      nParExt = nParExt / nParMesh                    
      IF( MODULO(nlev,nParExt) /= 0 ) THEN
        CALL Fatal(Caller,'Number of element layers '//I2S(nlev)//&
            ' not divisible by '//I2S(ParEnv % PEs))
      END IF
      nlev = nlev / nParExt
      IF(nlev < 2) THEN
        CALL Fatal(Caller,'At least two element layers needed in each partition!')
      END IF
      ilev = ( ParEnv % MyPe / nParMesh ) * nlev
      Wtable(0:nlev) = Wtable(ilev:nlev+ilev) 
    ELSE
      nParExt = 1
      nParMesh = 1 
      ilev = 0
    END IF
        
    ! Allocate extruded mesh:
    ! We do this only after splitting the division.
    ! ---------------------------------------------
    n = Mesh_in % NumberOfNodes
    nnodes = (nlev+1)*n

    Mesh_out => AllocateMesh()
    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )
    
    gnodes = Mesh_in % NumberOfNodes
    gelements = Mesh_in % NumberOfBulkElements

    Mesh_out % SingleMesh = .FALSE.
    SeparateSlices = .FALSE.
    CreateInternalBCs = .FALSE.
    
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

    CollectExtrudedBCs = ListGetLogical( CurrentModel % Simulation,'Extruded BCs Collect',Found )
    
    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo

      IF(.NOT. ASSOCIATED( PI_in ) ) CALL Fatal(Caller,'PI_in not associated!')
      IF(.NOT. ASSOCIATED( PI_out ) ) CALL Fatal(Caller,'PI_out not associated!')
            
      SeparateSlices = ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Slices Separate',Found )
      CreateInternalBCs = ListGetLogical( CurrentModel % Simulation,'Extruded BCs Internal',Found ) 
      IF(.NOT. Found) CreateInternalBCs = SeparateSlices 
      
      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % GInterface(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      IF(.NOT. SingleIn ) THEN
        IF(.NOT. ASSOCIATED( PI_in % NeighbourList ) ) THEN
          CALL Fatal(Caller,'Neighnours not associated in initial mesh!')
        END IF

        ! Count own nodes
        j=0
        DO i=1,Mesh_in % NumberOfNodes
          IF(.NOT. ASSOCIATED(PI_in % NeighbourList(i) % Neighbours ) ) THEN
            j = j + 1
          ELSE IF (PI_in % NeighbourList(i) % Neighbours(1) == ParEnv % MyPE ) THEN
            j=j+1
          END IF
        END DO
        CALL MPI_ALLREDUCE(j,gnodes,1, &
            MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
        gnodes = gnodes / nParExt
        
        j=0
        DO i=1,Mesh_in % NumberOfBulkElements
          IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
        END DO
        CALL MPI_ALLREDUCE(j,gelements,1, &
            MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
        gelements = gelements / nParExt

        !PRINT *,'nParExt:',ParEnv % Mype, nParExt, nParMesh, gnodes,gelements        
      END IF

      Mesh_out % ParallelInfo % NothingShared = ( SingleIn .AND. SeparateSlices )
    END IF

    CALL Info(Caller,'Number of nodes in layer: '//I2S(gnodes),Level=12)
    CALL Info(Caller,'Number of elements in layer: '//I2S(gelements),Level=12)
    
    cnt=0
    DO i=0,nlev

      CurrCoord = Wtable(i)  
      
      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          m = 1
          IF( nParMesh > 1 ) THEN
            IF( ASSOCIATED( PI_in % NeighbourList(j) % Neighbours ) ) THEN
              m = SIZE(PI_in % NeighbourList(j) % Neighbours)
            END IF
          END IF
          IF( SeparateSlices ) THEN
            k = m
          ELSE IF(i==0 .AND. ParEnv % MyPe > (nParMesh-1) ) THEN            
            k = 2*m
          ELSE IF(i==nlev .AND. ParEnv % MyPe < ParEnv % PEs- nParMesh ) THEN
            k = 2*m
          ELSE
            k = m
          END IF

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(k))
          PI_out % GInterface(cnt) = (k>1)
                    
          DO k=1,m
            IF(m>1) THEN
              OrigPart = PI_in % NeighbourList(j) % Neighbours(k)
            ELSE
              OrigPart = ParEnv % MyPe
            END IF                       

            IF( SeparateSlices ) THEN
              l = ( ParEnv % MyPe / nParMesh ) * ( nlev + 1 ) * gnodes
            ELSE
              l = ilev * gnodes
            END IF
                            
            IF(SingleIn) THEN
              l = l + j + i * gnodes
            ELSE
              l = l + MODULO(PI_in % GlobalDOFs(j)-1,gnodes)+1 + i * gnodes 
            END IF
            PI_out % GlobalDOFs(cnt) = l

            IF( SeparateSlices ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(k) = OrigPart               
            ELSE IF(i==0 .AND. ParEnv % MyPe > nParMesh-1 ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(2*k-1) = OrigPart
              PI_out % NeighbourList(cnt) % Neighbours(2*k) = OrigPart-1            
            ELSE IF(i==nlev .AND. ParEnv % MyPe < ParEnv % PEs-nParMesh ) THEN
              PI_out % NeighbourList(cnt) % Neighbours(2*k-1) = OrigPart+1
              PI_out % NeighbourList(cnt) % Neighbours(2*k) = OrigPart
            ELSE
              PI_out % NeighbourList(cnt) % Neighbours(k) = OrigPart 
            END IF                       
          END DO
          
        END IF
      END DO
    END DO
    
    Mesh_out % NumberOfNodes = cnt
    Mesh_out % Nodes % NumberOfNodes = cnt

    ! Calculate exactly and allocate the number of extruded elements
    n = Mesh_in % NumberOfBulkElements + Mesh_in % NumberOfBoundaryElements
    totalnumberofelements = n*nlev

    IF( CreateInternalBCs ) THEN
      totalnumberofelements = &
          totalnumberofelements + 2 * Mesh_in % NumberOfBulkElements 
    ELSE
      IF( ParEnv % MyPe < nParMesh ) totalnumberofelements = &
          totalnumberofelements + Mesh_in % NumberOfBulkElements 
      IF( ParEnv % MyPe >= ParEnv % PEs-nParMesh ) totalnumberofelements = &
          totalnumberofelements + Mesh_in % NumberOfBulkElements 
    END IF
      
    ALLOCATE(Mesh_out % Elements(totalnumberofelements))

    
    ! Generate volume bulk elements:
    ! ------------------------------
    Mesh_out % MaxElementNodes = 0
    n = Mesh_in % NumberOfNodes
    cnt=0; dg_n  = 0

    DO i=0,nlev-1
      DO j=1,Mesh_in % NumberOfBulkElements

        cnt = cnt+1
        Element => Mesh_out % Elements(cnt)
        Element = Mesh_in % Elements(j)

        bodyid = Element % BodyId
        Element % BodyId = MapExtrudedMaterial(Vlist,bodyid,ilev+i+1)        
        
        l_n = Mesh_in % Elements(j) % TYPE % NumberOfNodes
        ind(1:l_n) = Mesh_in % Elements(j) % NodeIndexes(1:l_n)+i*n
        ind(l_n+1:2*l_n) = Mesh_in % Elements(j) % NodeIndexes(1:l_n)+(i+1)*n
        l_n = 2*l_n
        Element % NDOFs = l_n
        Mesh_out % MaxElementNodes = MAX(Mesh_out % MaxElementNodes,l_n)

        SELECT CASE(l_n)
        CASE(6)
          Element % TYPE => GetElementType(706)
        CASE(8)
          Element % TYPE => GetElementType(808)
        END SELECT

        IF( isParallel ) THEN
          IF(SingleIn) THEN
            l = j + (ilev+i) * gelements
          ELSE
            l = MODULO(Mesh_in % Elements(j) % GElementIndex-1,gelements)+1 + (ilev+i) * gelements 
          END IF
          Element % GElementIndex = l
        ELSE
          Element % GElementIndex = cnt
        END IF
          
        Element % ElementIndex = cnt
        ALLOCATE(Element % NodeIndexes(l_n)) 
        Element % NodeIndexes = ind(1:l_n)
      END DO
    END DO
    Mesh_out % NumberOfBulkElements = cnt

    
    ! Add side boundaries with the bottom mesh boundary id's:
    ! -------------------------------------------------------
    max_bid = 0
    DO i=0,nlev-1
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Element => Mesh_out % Elements(cnt)
        Element = Mesh_in % Elements(k)        
        ALLOCATE(Element % BoundaryInfo)

        Element % BoundaryInfo = Mesh_in % Elements(k) % BoundaryInfo

        bid = Mesh_in % Elements(k) % BoundaryInfo % Constraint
        max_bid = MAX(max_bid, bid )

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l = Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Element % BoundaryInfo % Left => &
              Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l = Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Element % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        ElemCode = Mesh_in % Elements(k) % TYPE % ElementCode        
        m = 2*MODULO(ElemCode,100)        
        Element % NDOFs = m
        ALLOCATE(Element % NodeIndexes(m))
        pInds => Element % NodeIndexes
               
        IF(ElemCode == 202) THEN
          pInds(1) = Mesh_in % Elements(k) % NodeIndexes(1)+i*n
          pInds(2) = Mesh_in % Elements(k) % NodeIndexes(2)+i*n
          pInds(3) = Mesh_in % Elements(k) % NodeIndexes(2)+(i+1)*n
          pInds(4) = Mesh_in % Elements(k) % NodeIndexes(1)+(i+1)*n
          Mesh_out % Elements(cnt) % TYPE => GetElementType(404)
        ELSE IF(ElemCode == 101 ) THEN
          pInds(1) = Mesh_in % Elements(k) % NodeIndexes(1) +i*n
          pInds(2) = Mesh_in % Elements(k) % NodeIndexes(1) +(i+1)*n
        ELSE
          CALL Fatal(Caller,'Cannot extrude boundary element: '//I2S(ElemCode))
        END IF
        Element % ElementIndex = cnt
      END DO
    END DO

    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
   
    CALL Info(Caller,'First Extruded BC set to: '//I2S(max_bid+1),Level=6)
    lastbc = max_bid+1

    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    IF( CollectExtrudedBCs ) THEN
      CALL Info(Caller,'Number of new BCs for each layer: 1',Level=6)
    ELSE
      CALL Info(Caller,'Number of new BCs for each layer: '//I2S(max_body),Level=6)
    END IF
    
    IF( CollectExtrudedBCs ) THEN
      newbcs = 2
    ELSE
      newbcs = 2 * max_body
    END IF

    IF( CreateInternalBCs ) THEN
      CALL Info(Caller,'Internal bottom boundary: '//I2S(max_bid+newbcs+1),Level=6)
      CALL Info(Caller,'Internal top boundary: '//I2S(max_bid+newbcs+2),Level=6)
    END IF
    
    ALLOCATE(ChildBCs(2*max_body))
    ChildBCs = -1
           
    ! Add bottom boundary:
    ! --------------------
    IF( ParEnv % PEs == 1 .OR. ParEnv % MyPe < nParMesh .OR. CreateInternalBCs ) THEN  
      InternalBC = (ParEnv % PEs > 1 .AND. ParEnv % MyPe >= nParMesh )
      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        Element => Mesh_out % Elements(cnt) 
        
        Element = Mesh_in % Elements(i)

        l_n = Mesh_in % Elements(i) % TYPE % NumberOfNodes
        Element % NDOFs = l_n

        ALLOCATE(Element % BoundaryInfo)
        Element % BoundaryInfo % Left => Mesh_out % Elements(i)
        Element % BoundaryInfo % Right => NULL()

        bodyid = Mesh_in % Elements(i) % BodyId                
        IF( InternalBC ) THEN
          bcid = max_bid + newbcs + 1
        ELSE IF( CollectExtrudedBCs ) THEN
          bcid = max_bid + 1
        ELSE
          bcid = max_bid + bodyid
        END IF
        Element % BoundaryInfo % Constraint = bcid

        IF(.NOT. InternalBC) ChildBCs(2*bodyid-1) = bcid
        lastbc = MAX(lastbc,bcid)

        Element % BodyId = 0
        IF( bcid <= CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
          IF(Found) Element % BodyId = j
        END IF

        ALLOCATE(Element % NodeIndexes(l_n))
        Element % NodeIndexes = Mesh_in % Elements(i) % NodeIndexes
        Element % ElementIndex = cnt
        Element % TYPE => Mesh_in % Elements(i) % TYPE
      END DO
    END IF

    
    ! Add top boundary:
    ! -----------------
    IF( ParEnv % PEs == 1 .OR. ParEnv % MyPe >= ParEnv % PEs - nParMesh .OR. CreateInternalBCs ) THEN
      InternalBC = (ParEnv % PEs > 1 .AND. ParEnv % MyPe < ParEnv % PEs - nParMesh )
      DO i=1,Mesh_in % NumberOfBulkElements
        cnt=cnt+1
        Element => Mesh_out % Elements(cnt) 
        
        Element = Mesh_in % Elements(i)

        l_n = Mesh_in % Elements(i) % TYPE % NumberOfNodes
        Element % NDOFs = l_n

        ALLOCATE(Element % BoundaryInfo)
        Element % BoundaryInfo % Left => &
            Mesh_out % Elements((nlev-1)*Mesh_in % NumberOfBulkElements+i)
        Element % BoundaryInfo % Right => NULL()
        
        bodyid = Mesh_in % Elements(i) % BodyId                
        IF( InternalBC ) THEN
          bcid = max_bid + newbcs + 2
        ELSE IF( CollectExtrudedBCs ) THEN
          bcid = max_bid + 2
        ELSE
          bcid = max_bid + bodyid + max_body
        END IF
        Element % BoundaryInfo % Constraint = bcid

        IF(.NOT. InternalBC) ChildBCs(2*bodyid) = bcid 
        lastbc = MAX(lastbc,bcid)
        
        Element % BodyId = 0
        IF( bcid<=CurrentModel % NumberOfBCs) THEN
          j = ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
          IF(Found) Element % BodyId = j
        END IF

        ALLOCATE(Element % NodeIndexes(l_n))
        Element % NodeIndexes = Mesh_in % Elements(i) % NodeIndexes+nlev*n
        Element % ElementIndex = cnt
        Element % TYPE => Mesh_in % Elements(i) % TYPE
      END DO
    END IF

    IF(.NOT. SingleIn .AND. isParallel) THEN
      j=lastbc
      CALL MPI_ALLREDUCE(j,lastbc,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF
    CALL Info(Caller,'Last Extruded BC set to: '//I2S(lastbc),Level=6)
    
    IF( cnt /= totalnumberofelements ) THEN
      CALL Fatal(Caller,'Mismatch between allocated and set elements: '//&
          I2S(totalnumberofelements)//' vs. '//I2S(cnt))
    END IF

    ! Set some unset stuff to be on the safe side
    DO i=1,cnt
      Element => Mesh_out % Elements(i)
      Element % DGDOFs = 0
      Element % DGIndexes => NULL()
      Element % PDefs => NULL()
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()
      Element % BubbleIndexes => NULL()
    END DO
         
    Mesh_out % NumberOfBoundaryElements = cnt - Mesh_out % NumberOfBulkElements
    
    Mesh_out % Name = Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs = Mesh_out % MaxElementNodes
    Mesh_out % Stabilize = Mesh_in % Stabilize
    Mesh_out % MeshDim = 3
    CurrentModel % DIMENSION = 3


    ! Let us mark the child BCs to the bodies that they originate from.
    BLOCK
      INTEGER, POINTER :: TmpPair(:), TmpBCs(:) 
      TYPE(ValueList_t), POINTER :: vList

      ALLOCATE(TmpBCs(2*max_body))
      TmpBCs = ChildBCs

      IF( ParEnv % PEs > 1 ) THEN
        CALL MPI_ALLREDUCE(TmpBCs,ChildBCs,2*max_body, &
            MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
      END IF

      DO i=1,CurrentModel % NumberOfBodies
        vList => CurrentModel % Bodies(i) % Values
        IF( ASSOCIATED(vList) ) THEN
          NULLIFY(TmpPair)
          ALLOCATE(TmpPair(2))
          TmpPair(1) = ChildBCs(2*i-1)
          TmpPair(2) = ChildBCs(2*i)
          CALL ListAddIntegerArray(vList,'Extruded Child BCs',2,TmpPair)

          IF( InfoActive(10) ) THEN
            CALL Info(Caller,'Setting Body '//I2S(i)//' "Extruded Child BCs" to '&
                //I2S(TmpPair(1))//' '//I2S(TmpPair(2)))
          END IF
          NULLIFY(TmpPair)
        END IF
      END DO

      DEALLOCATE(TmpBCs)
    END BLOCK
      
    CALL SetMeshSkew(Mesh_out, CurrentModel % Simulation )
    
    ExtrudedMeshName = ListGetString(CurrentModel % Simulation,'Extruded Mesh Name',Found)
    IF(Found) THEN
      IF( ParEnv % PEs == 1 ) THEN
        CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
      ELSE
        CALL WriteMeshToDisk2(CurrentModel, Mesh_out, ExtrudedMeshName, ParEnv % MyPe )
      END IF
    END IF

    CALL PrepareMesh( CurrentModel, Mesh_out, isParallel )
    
!------------------------------------------------------------------------------
  END FUNCTION MeshExtrudeSlices
!------------------------------------------------------------------------------


  ! Routine for increasing element order by adding an additional node an each edge.
  ! Basically the same order of elements could be created by p-elements but this provides
  ! alternative solution when nodal finite element are preferred. Often the mesh may be
  ! made quadratic with the preprocessors but this enables also the use of mesh extrusion
  ! and mesh multiplication which cannot be used with higher order nodal elements.
  !--------------------------------------------------------------------------------------
  SUBROUTINE IncreaseElementOrder( Model, Mesh )
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
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
    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, NewMesh % NumberOfBoundaryElements
    
    WRITE( 1,'(i0)' ) 2
    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(k) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(k) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(k) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBoundaryElements

    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(i) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(i) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(i) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBulkElements
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.nodes', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       WRITE(1,'(i0,a,3e23.15)',ADVANCE='NO') i,' -1 ', &
            NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.elements', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       WRITE(1,'(3(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
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
       WRITE(1,'(5(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(k) % BoundaryInfo % Constraint, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
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
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: NewMesh
    CHARACTER(LEN=*) :: Path
    INTEGER, OPTIONAL :: Partition
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeList(100),ElmCodeCounts(100),&
        Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared, &
        iostat, BCWarns
    INTEGER, POINTER :: BList(:)
    INTEGER, ALLOCATABLE :: ElementCodes(:)
    LOGICAL :: Parallel, WarnNoTarget, Found
    CHARACTER(:), ALLOCATABLE :: headerFN, elementFN, nodeFN,&
         boundFN, sharedFN
!------------------------------------------------------------------------------

    IF(PRESENT(Partition)) THEN
       Parallel = .TRUE.
       headerFN = '/part.'//I2S(Partition+1)//'.header'
       elementFN = '/part.'//I2S(Partition+1)//'.elements'
       nodeFN =  '/part.'//I2S(Partition+1)//'.nodes'
       boundFN = '/part.'//I2S(Partition+1)//'.boundary'
       sharedFN ='/part.'//I2S(Partition+1)//'.shared'
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
    OPEN( 1,FILE=TRIM(Path) // headerFN,STATUS='UNKNOWN', iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//headerFN)
    END IF

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
    OPEN( 1,FILE=TRIM(Path) // nodeFN, STATUS='UNKNOWN',iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//nodeFN)
    END IF
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
    OPEN( 1,FILE=TRIM(Path) // elementFN, STATUS='UNKNOWN', iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//elementFN)
    END IF
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
    OPEN( 1,FILE=TRIM(Path) // boundFN, STATUS='UNKNOWN',iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//boundFN)
    END IF
    BcWarns = 0
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

       IF(.NOT. ASSOCIATED(NewMesh % Elements(k) % BoundaryInfo ) ) THEN
         CALL Fatal('WriteMeshToDisk2','BoundaryInfo not associated for element: '//I2S(k))
       END IF
       
       Constraint = NewMesh % Elements(k) % BoundaryInfo % Constraint

       Found = .FALSE.
       IF(Constraint > 0 .AND. Constraint <= Model % NumberOfBCs ) THEN
         BList => ListGetIntegerArray( Model % BCs(Constraint) % Values, &
             'Target Boundaries', Found )
       END IF
       IF(Found) THEN
          IF(SIZE(BList) > 1) THEN
            BcWarns = BcWarns + 1
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

    IF(BcWarns > 1 ) THEN
      CALL WARN("WriteMeshToDisk2",&
          "BC elements '//I2S(BcWarns)//' have more than one Target Boundary, SaveMesh output will not match input!")
    END IF
      
    IF(WarnNoTarget) THEN
       CALL WARN("WriteMeshToDisk2","Couldn't find a Target Boundary, assuming mapping to self")
    END IF

    IF(.NOT. Parallel) RETURN

    !Write .shared file
    !Need to create part.n.shared from Mesh % ParallelInfo %
    !NeighbourList % Neighbours.
    OPEN( 1,FILE=TRIM(Path) // sharedFN, STATUS='UNKNOWN',iostat=iostat)
    IF(iostat /= 0) THEN
      CALL Fatal('WriteMeshToDisk2','Could not open file: '//TRIM(Path)//sharedFN)
    END IF
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
    CHARACTER(:), ALLOCATABLE :: DirectoryName, PrefixName
!------------------------------------------------------------------------------

    NoPartitions = MAXVAL( ElementPart ) 
    NumElmCodes = 0
    NumElements = Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements
        
    DirectoryName = TRIM(PATH)//'/partitioning.'//I2S(NoPartitions)
    CALL MakeDirectory( DirectoryName // CHAR(0) )
    CALL Info('WriteMeshToDiskPartitioned','Writing parallel mesh to disk: '//DirectoryName)
   

    DO Partition = 1, NoPartitions 
      
      CALL Info('WriteMeshToDiskPartitioned','Writing piece to file: '//I2S(Partition),Level=12)
      
      PrefixName = DirectoryName//'/part.'//I2S(Partition)

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
!> Show mesh size information
!------------------------------------------------------------------------------
  SUBROUTINE PrintMeshSize( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: na, nb, nn, ne, nf, no, ns, i
    INTEGER :: napar(0:2), nbpar(0:2), nnpar(0:2), nepar(0:2), nfpar(0:2), nopar(0:2), nspar(0:2)
    CHARACTER(*), PARAMETER :: Caller="PrintMeshSize"   
!------------------------------------------------------------------------------

    na = Mesh % NumberOfBulkElements
    nb = Mesh % NumberOfBoundaryElements
    nn = Mesh % NumberOfNodes
    ne = Mesh % NumberOfEdges
    nf = Mesh % NumberOfFaces

    IF( ParEnv % PEs > 1 .AND. .NOT. Mesh % SingleMesh ) THEN
      no = 0; ns = 0
      DO i=1,nn
        IF(Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % MyPe) no = no+1
        IF(SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours) > 1) ns = ns+1
      END DO
      DO i=0,2
        napar(i) = ParallelReduction(na,i)
        nbpar(i) = ParallelReduction(nb,i)
        nnpar(i) = ParallelReduction(nn,i)
        nopar(i) = ParallelReduction(no,i)
        nspar(i) = ParallelReduction(ns,i)
        nepar(i) = ParallelReduction(ne,i)
        nfpar(i) = ParallelReduction(nf,i)
      END DO

      CALL Info(Caller,'Number of parallel mesh entities:   SUM       MIN       MAX')
      WRITE(Message,'(A,T30,3I10)') '  Bulk elements: ',napar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Boundary elements: ',nbpar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Total nodes: ',nnpar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Owned nodes: ',nopar
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,3I10)') '  Shared nodes: ',nspar
      CALL Info(Caller,Message,Level=3)
      IF(nepar(0) > 0) THEN
        WRITE(Message,'(A,T30,3I10)') '  Element edges: ',nepar
        CALL Info(Caller,Message,Level=3)
      END IF
      IF(nfpar(0) > 0) THEN
        WRITE(Message,'(A,T30,3I10)') '  Element faces: ',nfpar
        CALL Info(Caller,Message,Level=3)
      END IF
    ELSE
      CALL Info(Caller,'Number of serial mesh entities')
      WRITE(Message,'(A,T30,1I10)') '  Bulk elements: ',na
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,1I10)') '  Boundary elements: ',nb
      CALL Info(Caller,Message,Level=3)
      WRITE(Message,'(A,T30,1I10)') '  Element nodes: ',nn
      CALL Info(Caller,Message,Level=3)
      IF(ne > 0) THEN
        WRITE(Message,'(A,T30,1I10)') '  Element edges: ',ne
        CALL Info(Caller,Message,Level=3)
      END IF
      IF(nf > 0) THEN
        WRITE(Message,'(A,T30,1I10)') '  Element faces: ',nf
        CALL Info(Caller,Message,Level=3)
      END IF
    END IF

  END SUBROUTINE PrintMeshSize

      
  
!------------------------------------------------------------------------------
!> Check mesh for various info. Mainly for debugging.
!------------------------------------------------------------------------------
  SUBROUTINE CheckMeshInfo( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: na, nb, nn
    INTEGER :: i,j,k,t,ii,jj,maxi,mini
    INTEGER, ALLOCATABLE :: NodeHits(:), TypeHits(:)
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: mins, maxs, s, s2
    INTEGER(KIND=8) :: Dbg(10)
    LOGICAL :: Halt
    CHARACTER(*), PARAMETER :: Caller="CheckMeshInfo"   
!------------------------------------------------------------------------------

    CALL Info(Caller,'Checking mesh information')

    na = Mesh % NumberOfBulkElements
    nb = Mesh % NumberOfBoundaryElements
    nn = Mesh % NumberOfNodes
    Halt = .FALSE.
    
    ALLOCATE(TypeHits(827))
    ALLOCATE(NodeHits(nn))
    
    CALL CheckMeshBulkHits()
    CALL CheckMeshBoundaryHits()
    CALL CheckBCTags()
    CALL CheckParentIndeces()
    CALL CheckMeshGeomSize()
    CALL CheckMeshSerendipity()
    CALL CheckMeshBodyRadius()
    CALL CheckMeshFaces()
    CALL CheckMeshEdges()
    CALL CheckParallelInfo()
    CALL CheckParallelEdgeInfo()
    CALL CheckParallelFaceInfo()

    nn = ParallelReduction(nn)
    
    CALL Info(Caller,'Finished checking mesh!')

    IF(Halt) CALL Fatal(Caller,'Some checksum was invalid, cannot continue!')


  CONTAINS

    
    SUBROUTINE CheckMeshBulkHits()
      TypeHits = 0
      NodeHits = 0
      Dbg = 0

      DO t=1,na
        Element => Mesh % Elements(t)
        IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
          CALL Fatal(Caller,'Element type not associated for bulk elem: '//I2S(t))
        END IF
        i = Element % TYPE % ElementCode        
        TypeHits(i) = TypeHits(i)+1
        IF(ANY(Element % NodeIndexes < 1 ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes 
          CALL Fatal(Caller,'Bulk element '//I2S(t)//' has non-positive index!')
        END IF
        IF(ANY(Element % NodeIndexes > nn ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes, ' vs. ', nn 
          CALL Fatal(Caller,'Bulk element '//I2S(t)//' has too large index!')
        END IF
        IF(ANY(Element % NodeIndexes <= 0)) THEN
          PRINT *,'NodeIndexes:',Element % NodeIndexes
          CALL Fatal(Caller,'Non-positive node index encountered')
        END IF
        IF(ANY(Element % NodeIndexes < 1) ) THEN
          PRINT *,'Too small bulk element index: ',t, Element % NodeIndexes
        ELSE IF(ANY(Element % NodeIndexes > nn) ) THEN
          PRINT *,'Too large bulk element index: ',t, Element % NodeIndexes
        ELSE
          NodeHits(Element % NodeIndexes) = NodeHits(Element % NodeIndexes) + 1
        END IF
      END DO

      Dbg(1) = na
      Dbg(2) = SUM(NodeHits) 
      DO i=1,SIZE(NodeHits)
        Dbg(3) = dbg(3) + i*NodeHits(i)
      END DO

      DO i=1,SIZE(TypeHits)
        j = TypeHits(i)
        IF(j>0) CALL Info(Caller,'Bulk element type '//I2S(i)//' count: '//I2S(j))        
      END DO
      
      t=MAXVAL(NodeHits)

      IF( InfoActive(25)) THEN
        DO i=0,t
          j = COUNT( NodeHits == i)
          IF(j>0) PRINT *,'Bulk node hits '//I2S(i)//' count: ',j
        END DO
      END IF
      dbg(4) = t      
      Dbg(5) = COUNT(TypeHits>0)


      WRITE(Message,*) 'Bulk Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)      
      
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckMeshBulkHits


    SUBROUTINE CheckMeshBoundaryHits()
      TypeHits = 0
      NodeHits = 0
      dbg = 0

      DO t=na+1,na+nb
        Element => Mesh % Elements(t)
        IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
          CALL Fatal(Caller,'Element type not associated for bc elem: '//I2S(t-na))
        END IF
        i = Element % TYPE % ElementCode
        TypeHits(i) = TypeHits(i)+1
        IF(ANY(Element % NodeIndexes < 1 ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes 
          CALL Fatal(Caller,'Boundary element '//I2S(t)//' has non-positive index!')
        END IF
        IF(ANY(Element % NodeIndexes > nn ) ) THEN
          PRINT *,'NodeIndexes:', Element % NodeIndexes, ' vs. ', nn 
          CALL Fatal(Caller,'Boundary element '//I2S(t)//' has too large index!')
        END IF
        NodeHits(Element % NodeIndexes) = NodeHits(Element % NodeIndexes) + 1

        IF(ASSOCIATED(Element % BoundaryInfo) ) THEN
          IF(.NOT. ( ASSOCIATED(Element % BoundaryInfo % Left) .OR. &
              ASSOCIATED(Element % BoundaryInfo % Right) ) ) THEN
            PRINT *,'Boundary Info present but no left/right parent: ',t
          END IF
        END IF
      END DO
      DO i=1,SIZE(TypeHits)
        j = TypeHits(i)
        IF(j>0) CALL Info(Caller,'Boundary element type '//I2S(i)//' count: '//I2S(j))        
      END DO

      t=MAXVAL(NodeHits)
      IF(InfoActive(25)) THEN
        DO i=0,t
          j = COUNT( NodeHits == i)
          IF(j>0) PRINT *,'Boundary node hits '//I2S(i)//' count: ',j
        END DO
      END IF
        
      Dbg(1) = nb
      Dbg(2) = SUM(NodeHits) 
      DO i=1,SIZE(NodeHits)
        Dbg(3) = dbg(3) + i*NodeHits(i)
      END DO
      Dbg(4) = COUNT(TypeHits>0)
      Dbg(5) = t
      WRITE(Message,*) 'Boundary Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)      
     
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckMeshBoundaryHits


    SUBROUTINE CheckParentIndeces()
      INTEGER :: Misses
      TYPE(Element_t), POINTER :: Parent

      Misses = 0
      dbg = 0
      
      DO t=na+1,na+nb
        Element => Mesh % Elements(t)
        i = Element % TYPE % NumberOfNodes
        IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE
        DO j=1,2
          IF(j==1) THEN
            Parent => Element % BoundaryInfo % Left
          ELSE
            Parent => Element % BoundaryInfo % Right
          END IF
          IF(.NOT. ASSOCIATED(Parent)) CYCLE

          dbg(3) = dbg(3) + 1
          dbg(4) = dbg(4) + Element % ElementIndex
          dbg(5) = dbg(5) + Element % BoundaryInfo % Constraint

          k = 0
          DO i=1,Element % TYPE % NumberOfNodes
            IF( .NOT. ANY(Parent % NodeIndexes == Element % NodeIndexes(i) ) ) k=k+1
          END DO
          IF( k > 0 ) THEN
            Misses = Misses + 1
            IF( Misses <= 10 ) THEN
              PRINT *,'Indeces missing in parent: ',ParEnv % Mype, Element % ElementIndex,Parent % ElementIndex, &
                  Element % TYPE % NumberOfNodes, k
              PRINT *,'Element codes:',Element % TYPE % ElementCode, &
                  Parent % TYPE % elementCode
              PRINT *,'bc elem inds:',Element % NodeIndexes
              PRINT *,'bulk elem inds:',Parent % NodeIndexes 
            END IF
          END IF
        END DO
      END DO

      IF(Misses>0) PRINT *,'Parent elements missing nodes:',ParEnv % Mype, Misses      
      dbg(1) = nb
      dbg(2) = Misses
      
      WRITE(Message,*) 'Parent Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)

      IF(Misses > 0) CALL Fatal(Caller,'We need all parent indeces!')

      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckParentIndeces


    SUBROUTINE CheckBCTags()
      INTEGER :: Misses, Tag, MinTag, MaxTag
      INTEGER, ALLOCATABLE :: TagCount(:), BCNodeCount(:)

      ! Not BCs to go through
      IF(nb==0) RETURN

      Misses = 0
      MinTag = HUGE(MinTag)
      MaxTag = -HUGE(MaxTag)

      DO k=1,2
        DO t=na+1,na+nb
          Element => Mesh % Elements(t)
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) THEN
            IF(k==1) Misses = Misses + 1
            CYCLE
          END IF
          Tag = Element % BoundaryInfo % Constraint
          IF(k==1) THEN
            MinTag = MIN(MinTag,Tag)
            MaxTag = MAX(MaxTag,Tag)
          ELSE
            TagCount(Tag) = TagCount(Tag) + 1
          END IF
        END DO
        IF(k==1) THEN
          ! Not tags defined in this partition
          IF(MinTag > MaxTag) THEN
            EXIT
          ELSE
            ALLOCATE(TagCount(MinTag:MaxTag))
            TagCount = 0
          END IF
        END IF
      END DO

      ALLOCATE(BCNodeCount(MinTag:MaxTag))
      BCNodeCount = 0
      
      DO k=1, MaxTag
        IF(TagCount(k)==0) CYCLE
        NodeHits = 0
        DO t=na+1,na+nb
          Element => Mesh % Elements(t)
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE
          Tag = Element % BoundaryInfo % Constraint
          IF(Tag==k) NodeHits(Element % NodeIndexes) = 1
        END DO
        BCNodeCount(k) = COUNT(NodeHits == 1)
      END DO

      DO k=MinTag,MaxTag
        IF(TagCount(k) > 0) THEN
          PRINT *,'BC'//I2S(k)//': elems '//I2S(TagCount(k))//' nodes '//I2S(BCNodeCount(k))
        END IF
      END DO           
      
    END SUBROUTINE CheckBCTags

    
    
    SUBROUTINE CheckMeshGeomSize()

      IF(.NOT. InfoActive(25)) RETURN
      
      PRINT *,'Coordinate x: ',MINVAL(Mesh % Nodes % x), MAXVAL(Mesh % Nodes % x)
      PRINT *,'Coordinate y: ',MINVAL(Mesh % Nodes % y), MAXVAL(Mesh % Nodes % y)
      PRINT *,'Coordinate z: ',MINVAL(Mesh % Nodes % z), MAXVAL(Mesh % Nodes % z)

      mins = HUGE(mins); maxs = 0.0_dp
      DO t=1,na+nb
        Element => Mesh % Elements(t)
        DO i=1,Element % TYPE % NumberOfNodes
          ii = Element % NodeIndexes(i)
          DO j=i+1, Element % TYPE % NumberOfNodes
            jj = Element % NodeIndexes(j)
            s2 = (Mesh % Nodes % x(ii)-Mesh % Nodes % x(jj))**2 + &
                (Mesh % Nodes % y(ii)-Mesh % Nodes % y(jj))**2 + &
                (Mesh % Nodes % z(ii)-Mesh % Nodes % z(jj))**2
            IF( s2 < mins ) THEN
              mins = s2
              mini = t 
            END IF
            IF( s2 > maxs ) THEN
              maxs = s2
              maxi = t
            END IF
          END DO
        END DO

        IF( t==na .OR. t==na+nb) THEN
          mins = SQRT(mins)
          maxs = SQRT(maxs)            
          IF(t==na) THEN
            PRINT *,'Bulk element h range:',mins,maxs          
            mins = HUGE(mins); maxs = 0.0_dp
          ELSE
            PRINT *,'Boundary element h range:',mins,maxs
          END IF

          Element => Mesh % Elements(maxi)
          PRINT *,'Maximum element:',maxi
          PRINT *,'x:',Mesh % Nodes % x(Element % NodeIndexes)
          PRINT *,'y:',Mesh % Nodes % y(Element % NodeIndexes)
          PRINT *,'z:',Mesh % Nodes % z(Element % NodeIndexes)

        END IF
      END DO
      
    END SUBROUTINE CheckMeshGeomSize

    SUBROUTINE CheckMeshSerendipity()
      INTEGER :: ElemCode 
      INTEGER :: Indexes0(27),EdgeInds(2),n,ne
      INTEGER, POINTER :: Indexes(:)
      REAL(KIND=dp) :: Coord(3),Coord0(3)
      
      DO t=1,na
        Element => Mesh % Elements(t)

        n = Element % Type % NumberOfNodes
        ne = Element % Type % NumberOfEdges
        
        ElemCode = Element % TYPE % ElementCode
        Indexes => Element % NodeIndexes
        Indexes0(1:n) = Indexes(1:n)
        
        SELECT CASE( ElemCode )
        CASE( 306, 408 )

          DO i=1,ne
            EdgeInds(1) = Indexes(i)
            IF(i==ne) THEN
              EdgeInds(2) = Indexes(1)
            ELSE
              EdgeInds(2) = Indexes(i+1)
            END IF

            ! Center of edge 
            Coord0(1) = SUM( Mesh % Nodes % x(EdgeInds)) / 2
            Coord0(2) = SUM( Mesh % Nodes % y(EdgeInds)) / 2
            Coord0(3) = SUM( Mesh % Nodes % z(EdgeInds)) / 2

            ! Is there some node closer to center of edge?
            maxs = HUGE(maxs)
            DO j=ne+1,n
              Coord(1) = Mesh % Nodes % x(Indexes(j)) 
              Coord(2) = Mesh % Nodes % y(Indexes(j)) 
              Coord(3) = Mesh % Nodes % z(Indexes(j)) 
              s2 = SUM((Coord-Coord0)**2)
              IF(s2 < maxs ) THEN
                Indexes0(ne+i) = Indexes(j)
                maxs = s2
              END IF              
            END DO
          END DO

        END SELECT
          
        j = COUNT( Indexes(1:n) /= Indexes0(1:n) )
        IF( j > 0 ) THEN
          !PRINT *,'Discrepancy: ',Indexes(ne+1:n), Indexes0(ne+1:n)
          Element % NodeIndexes(1:n) = Indexes0(1:n)
          CALL Warn('CheckMeshInfo','Node order wrong for '//I2S(j)//' nodes in element '//I2S(t))
        END IF
          
      END DO
    END SUBROUTINE CheckMeshSerendipity


    SUBROUTINE CheckMeshBodyRadius()
      REAL(KIND=dp) :: r
      REAL(KIND=dp), ALLOCATABLE :: RadRange(:,:)
      INTEGER, POINTER :: Indexes(:)
      INTEGER, ALLOCATABLE :: BodyHits(:)
      INTEGER :: nob

      IF(.NOT. InfoActive(25)) RETURN
      
      nob = CurrentModel % NumberOfBodies
      ALLOCATE(RadRange(0:nob,2),BodyHits(0:nob))
      RadRange(:,1) = HUGE(r)
      RadRange(:,2) = 0.0_dp
      BodyHits = 0
      
      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        Indexes => Element % NodeIndexes
        r = 0.0_dp
        k = Element % BodyId
        IF(k<0 .OR. k>nob) CYCLE
        DO i=1, Element % TYPE % NumberOfNodes
          j = Indexes(i)
          r = Mesh % Nodes % x(j)**2
          r = r + Mesh % Nodes % y(j)**2
          r = r + Mesh % Nodes % z(j)**2
        END DO
        RadRange(k,1) = MIN(RadRange(k,1),r)
        RadRange(k,2) = MAX(RadRange(k,2),r)                   
        BodyHits(k) = BodyHits(k)+1
      END DO
      RadRange = SQRT( RadRange )
      DO i=0,nob
        IF(BodyHits(i)==0) CYCLE
        PRINT *,'Radius range: ',i,RadRange(i,:)
      END DO
    END SUBROUTINE CheckMeshBodyRadius


    SUBROUTINE CheckMeshEdges()
      INTEGER, POINTER :: Indexes(:)
      INTEGER :: m

      IF(Mesh % NumberOfEdges == 0 ) RETURN      
      dbg = 0
      dbg(1) = Mesh % NumberOfEdges

      DO t=1,Mesh % NumberOfEdges
        Element => Mesh % Edges(t)
        IF(.NOT. ASSOCIATED(Element)) THEN
          CALL Fatal(Caller,'Edge not associated on edge list: '//I2S(t))
        END IF
        Indexes => Element % NodeIndexes          
        IF(.NOT. ASSOCIATED(Indexes)) THEN
          CALL Fatal(Caller,'NodeIndexes not associated on edge: '//I2S(t))
        END IF
        IF(.NOT. ASSOCIATED(Element % TYPE)) THEN
          CALL Fatal(Caller,'Edge type '//I2S(t)//' not associated!')
        END IF
        m = Element % Type % NumberOfNodes
        IF(SIZE(Indexes) /= m) THEN
          CALL Fatal(Caller,'Invalid size of edge '//I2S(t)//&
              ' NodeIndexes: '//I2S(SIZE(Indexes))//' vs. '//I2S(m))
        END IF
        IF(SIZE(Indexes)>0) dbg(2) = dbg(2) + SUM(Indexes)
        dbg(3) = dbg(3) + Element % ElementIndex
        dbg(4) = dbg(4) + Element % GElementIndex        
      END DO

      WRITE(Message,*) 'Edges Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)                  

      !IF(ANY(Dbg < 0) ) Halt = .TRUE.

    END SUBROUTINE CheckMeshEdges

    SUBROUTINE CheckMeshFaces()
      INTEGER, POINTER :: Indexes(:)
      INTEGER :: m

      IF(Mesh % NumberOfFaces == 0 ) RETURN      
      dbg = 0
      dbg(1) = Mesh % NumberOfFaces

      DO t=1,Mesh % NumberOfFaces
        Element => Mesh % Faces(t)
        IF(.NOT. ASSOCIATED(Element)) THEN
          CALL Fatal(Caller,'Face not associated on face list: '//I2S(t))
        END IF
        Indexes => Element % NodeIndexes          
        IF(.NOT. ASSOCIATED(Indexes)) THEN
          CALL Fatal(Caller,'NodeIndexes not associated on face: '//I2S(t))
        END IF
        IF(.NOT. ASSOCIATED(Element % TYPE)) THEN
          CALL Fatal(Caller,'Face type '//I2S(t)//' not associated!')
        END IF
        m = Element % TYPE % NumberOfNodes
        IF(SIZE(Indexes) /= m) THEN
          CALL Fatal(Caller,'Invalid size of face '//I2S(t)//&
              ' NodeIndexes: '//I2S(SIZE(Indexes))//' vs. '//I2S(m))
        END IF
        IF(SIZE(Indexes)>0) dbg(2) = dbg(2) + SUM(Indexes)
        dbg(3) = dbg(3) + Element % ElementIndex
        dbg(4) = dbg(4) + Element % GElementIndex        

      END DO

      WRITE(Message,*) 'Faces Checksum: ',Dbg(1:5)
      CALL Info(Caller,Message)                  

      !IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckMeshFaces


    SUBROUTINE CheckParallelInfo()

      IF( ParEnv % PEs == 1) RETURN
      IF(.NOT. ASSOCIATED( Mesh % ParallelInfo % NeighbourList) ) RETURN
      
      dbg = 0
      dbg(1) = SIZE(Mesh % ParallelInfo % NeighbourList)

      dbg(2) = COUNT(Mesh % ParallelInfo % Ginterface)
      DO i=1, SIZE(Mesh % ParallelInfo % Ginterface)        
        IF( Mesh % ParallelInfo % Ginterface(i) ) dbg(3) = dbg(3) + i 
      END DO
      
      DO i=1, SIZE(Mesh % ParallelInfo % NeighbourList)
        IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) THEN
          dbg(7) = dbg(7) + 1
          CYCLE
        END IF
        j = SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
        dbg(4) = dbg(4) + j
        dbg(5) = dbg(5) + SUM(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
        dbg(6) = dbg(6) + i*j
      END DO

      WRITE(Message,*) 'ParallelInfo Checksum: ',Dbg(1:7)
      CALL Info(Caller,Message)                         

      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
     END SUBROUTINE CheckParallelInfo
       
    SUBROUTINE CheckParallelEdgeInfo()

      IF( ParEnv % PEs == 1) RETURN
      IF( Mesh % NumberOfEdges == 0) RETURN
      IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % EdgeNeighbourList)) RETURN
      
      dbg = 0      
      dbg(1) = SIZE(Mesh % ParallelInfo % EdgeNeighbourList)

      IF(ASSOCIATED(Mesh % ParallelInfo % EdgeInterface ) ) THEN
        j = SIZE(Mesh % ParallelInfo % EdgeInterface )        
        IF(j>1) THEN
          dbg(2) = j
          DO i=1, SIZE(Mesh % ParallelInfo % Edgeinterface)        
            IF( Mesh % ParallelInfo % Edgeinterface(i) ) dbg(3) = dbg(3) + i 
          END DO
        END IF
      END IF
        
      DO i=1, SIZE(Mesh % ParallelInfo % EdgeNeighbourList)
        IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)) THEN
          dbg(7) = dbg(7) + 1
          CYCLE
        END IF
        j = SIZE(Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)
        dbg(4) = dbg(4) + j
        dbg(5) = dbg(5) + SUM(Mesh % ParallelInfo % EdgeNeighbourList(i) % Neighbours)
        dbg(6) = dbg(6) + i*j
      END DO

      WRITE(Message,*) 'ParallelEdges Checksum: ',Dbg(1:7)
      CALL Info(Caller,Message)                         
      
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckParallelEdgeInfo

    SUBROUTINE CheckParallelFaceInfo()

      IF( ParEnv % PEs == 1) RETURN
      IF( Mesh % NumberOfFaces == 0) RETURN
      IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % FaceNeighbourList)) RETURN
      
      dbg = 0
      dbg(1) = SIZE(Mesh % ParallelInfo % FaceNeighbourList)

      IF( ASSOCIATED( Mesh % ParallelInfo % FaceInterface ) ) THEN
        j = SIZE(Mesh % ParallelInfo % FaceInterface )
        IF(j>1) THEN
          dbg(2) = j
          DO i=1, j
            IF( Mesh % ParallelInfo % Faceinterface(i) ) dbg(3) = dbg(3) + i 
          END DO
        END IF
      END IF

      DO i=1, SIZE(Mesh % ParallelInfo % FaceNeighbourList)
        IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)) THEN
          dbg(7) = dbg(7) + 1
          CYCLE
        END IF
        j = SIZE(Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)
        dbg(4) = dbg(4) + j
        dbg(5) = dbg(5) + SUM(Mesh % ParallelInfo % FaceNeighbourList(i) % Neighbours)
        dbg(6) = dbg(6) + i*j
      END DO

      WRITE(Message,*) 'ParallelFaces Checksum: ',Dbg(1:7)
      CALL Info(Caller,Message)                         
      
      IF(ANY(Dbg < 0) ) Halt = .TRUE.
      
    END SUBROUTINE CheckParallelFaceInfo
           
!------------------------------------------------------------------------------
  END SUBROUTINE CheckMeshInfo
!------------------------------------------------------------------------------


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
    TYPE(Mesh_t), POINTER :: Mesh
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
    TYPE(Mesh_t), POINTER :: Mesh
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
    LOGICAL :: Parallel
    CHARACTER(*), PARAMETER :: Caller = 'SplitMeshEqual'
!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Mesh ) ) RETURN

    CALL Info( Caller, 'Mesh splitting works for first order elements 303, 404, 504, (706) and 808.', Level = 6 )

    DO i=1,Mesh % NumberOfBulkElements
      SELECT CASE(Mesh % Elements(i) % TYPE % ElementCode/100)
      CASE(6)
        CALL Fatal(Caller,'Pyramids not supported, sorry.')
      END SELECT
    END DO

    NewMesh => AllocateMesh()

    NewMesh % SingleMesh = Mesh % SingleMesh
    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. NewMesh % SingleMesh )

    
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT.EdgesPresent) CALL FindMeshEdges( Mesh )

    CALL ResetTimer(Caller)

    CALL Info( Caller, '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( Caller, Message, Level=6 )
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
    IF(FaceCnt>0) CALL Info( Caller,'Added '//I2S(FaceCnt)//' nodes in the center of faces',Level=10)

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
    IF(NodeIt>0) CALL Info( Caller,'Added '//I2S(NodeIt)//' nodes in the center of bulks',Level=10)

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
    x(1:Mesh % NumberOfNodes) = u(1:Mesh % NumberOfNodes)
    y(1:Mesh % NumberOfNodes) = v(1:Mesh % NumberOfNodes)
    z(1:Mesh % NumberOfNodes) = w(1:Mesh % NumberOfNodes)

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
    CALL Info(Caller,'Added edge centers to the nodes list.', Level=15 )  

!   add quad face centers for bricks and prisms(wedges):
!   ----------------------------
    j = Mesh % NumberOfNodes + Mesh % NumberOfEdges
    DO i=1,Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       k = Face % TYPE % NumberOfNodes
       IF( k==4 ) THEN
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
    CALL Info(Caller,'Added face centers to the nodes list.', Level=15 )

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
    CALL Info(Caller,'Added quad and brick centers to the nodes list.', Level=15 )

    
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MinEdgeDOFs = Mesh % MinEdgeDOFs
    NewMesh % MinFaceDOFs = Mesh % MinFaceDOFs
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
       CASE(1)
          NewElCnt = NewElCnt + 1 ! lines
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
    CALL Info( Caller, Message, Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info(Caller,'New mesh allocated.', Level=10 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 8 )
    CALL Info(Caller,'Array for bulk elements allocated.', Level=10 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes
    EdgeCnt = Mesh % NumberOfEdges

!
!   Index to old edge/quad/hexa centerpoint node in the new mesh nodal arrays:
!   ---------------------------------------------------------------------
    Node = NodeCnt + EdgeCnt + FaceCnt
!
!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)

       SELECT CASE( Eold % TYPE % ElementCode )
       CASE(101)
!
!         Copy point element
!         ------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 1)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)

       CASE(202)
!
!         Split edge to two edges
!         ------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 2)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,2) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 2)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)

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
              ' not supported by the multigrid solver.'
          CALL Fatal( Caller, Message )
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
       CASE(1)
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 1 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
          DO j=NewElCnt-1,1,-1
            Eold => NewMesh % Elements(j)
            IF(Eold % Type % ElementCode/100==2) THEN
              IF(ANY(Eold % NodeIndexes==Enew % NodeIndexes(1))) THEN
                IF(.NOT.ASSOCIATED(Enew % BoundaryInfo % Left)) THEN
                  Enew % BoundaryInfo % Left => Eold
                ELSE
                  Enew % BoundaryInfo % Right => Eold
                  EXIT
                END IF
              END IF
            END IF
          END DO

       CASE(2)
!
!         Line segments:
!         ==============
!
!         which edge of the parent element are we ?
!         -----------------------------------------
          Found = .FALSE.
          DO Edge1=1,SIZE(Eparent % EdgeIndexes)
            Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
            Found = ANY(Eold % NodeIndexes(1:2) == Edge % NodeIndexes(1) ) .AND. &
                ANY(Eold % NodeIndexes(1:2) == Edge % NodeIndexes(2) )
            IF(Found) EXIT
          END DO
          IF(.NOT. Found) THEN
            CALL Fatal(Caller,'Could not find parent edge with nodes: '//&
                I2S(Eold % NodeIndexes(1))//' '//I2S(Eold % NodeIndexes(2)))
          END IF

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

          Found = .FALSE.

          n1 = 4 
          IF( Eparent % TYPE % ElementCode > 500 ) n1 = 8
          
          DO j=1,n1
            Eptr => NewMesh % Elements( Child(ParentId,j) )
            n = Eptr % TYPE % NumberOfNodes

            ! The parent is unique! Hence it is enough to find a parent with both matches.
            Found =  ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(1) ) .AND. &
                ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(2) )
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
                    
          DO j=1,n1
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found =  ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(1) ) .AND. &
                 ANY( Eptr % NodeIndexes(1:n) == Enew % NodeIndexes(2) )
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
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
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
               IF( ANY( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
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
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
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
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
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
               IF( ANY( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n) ) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 )  CALL Error( Caller, 'Parent element not found' )
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
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
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
               IF( ANY( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
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
               IF( ANY(Enew % NodeIndexes(n1) == Eptr % NodeIndexes(1:n)) ) n3 = n3+1
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( Caller, 'Parent element not found' )
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
           Enew % PDefs % isEdge = .FALSE.

           ! If element is of type tetrahedron and is a p element,
           ! do the Ainsworth & Coyle trick
           IF (Enew % TYPE % ElementCode == 504) CALL ConvertToACTetra(Enew)
            CALL GetRefPElementNodes( Enew % Type,  Enew % Type % NodeU, &
                 Enew % Type % NodeV, Enew % Type % NodeW )
         END IF
      ELSE
        Enew % PDefs=>NULL()
      END IF
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()
      Enew % BubbleIndexes => NULL()
    END DO

    CALL Info( Caller, '******** New mesh ********', Level=6 )
    WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
    CALL Info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
    CALL Info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
    CALL Info( Caller, Message, Level=6 )


    ! Information of the new system size, also in parallel
    !----------------------------------------------------------------------
    ParTmp(1) = Mesh % NumberOfNodes
    ParTmp(2) = Mesh % NumberOfBulkElements
    ParTmp(3) = Mesh % NumberOfBoundaryElements
    ParTmp(4) = NewMesh % NumberOfNodes
    ParTmp(5) = NewMesh % NumberOfBulkElements
    ParTmp(6) = NewMesh % NumberOfBoundaryElements

    IF( .FALSE. .AND. Parallel ) THEN
      CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)

      CALL Info(Caller,'Information on parallel mesh sizes')
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(1),' nodes'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(2),' bulk elements'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(3),' boundary elements'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(4),' nodes'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(5),' bulk elements'
      CALL Info(Caller,Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(6),' boundary elements'
      CALL Info(Caller,Message)
    END IF


    CALL CheckTimer(Caller,Delete=.TRUE.)

!
!   Update structures needed for parallel execution:
!   ------------------------------------------------
    IF( Parallel ) THEN
      CALL UpdateParallelMesh( Mesh, NewMesh )
    END IF
!
!
!   Finalize:
!   ---------
    DEALLOCATE( Child )
    IF(.NOT.EdgesPresent) THEN
      CALL Info(Caller,'Releasing edges from the old mesh as they are not needed!',Level=20)
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    ELSE
      CALL Info(Caller,'Generating edges in the new mesh as they were present in the old!',Level=20)
      CALL FindMeshEdges( NewMesh )
    END IF

    ! Our boundary may be a circle, cylinder or sphere surface.
    ! Honor those shapes when splitting the mesh!
    CALL FollowCurvedBoundary( CurrentModel, NewMesh, .FALSE. ) 
    
    
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
       CALL AllocateVector( NewMesh % ParallelInfo % GInterface,n  )
       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )

       DO i=1,n
          NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % GInterface = .FALSE.
       NewMesh % ParallelInfo % GInterface(1:n) = Mesh % ParallelInfo % GInterface

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
             IF( Mesh % ParallelInfo % GInterface(i) ) p = p+1
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
             IF( .NOT.ALL( Mesh % ParallelInfo % GInterface( Edge % NodeIndexes ) )) CYCLE
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
             IF( NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + i) ) CYCLE
             
             ! Mark interface nodes and count edges:
             !--------------------------------------
             NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + i) = .TRUE.
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
             IF( Mesh % ParallelInfo % GInterface(i) ) p = p+1
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
             IF( .NOT.ALL( Mesh % ParallelInfo % GInterface( Face % NodeIndexes ) )) CYCLE
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
                IF( NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + k) ) CYCLE
                
                ! Mark interface nodes and count edges:
                !--------------------------------------
                NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes + k) = .TRUE.
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
             IF ( ALL( Mesh % ParallelInfo % GInterface( Face % NodeIndexes ) ) ) THEN
                NewMesh % ParallelInfo % GInterface( Mesh % NumberOfNodes &
                     + Mesh % NumberOfEdges + i ) = .TRUE.
                j = j + 1
                k = k + Face % TYPE % NumberOfNodes
             END IF
          END IF
       END DO

!      CALL AllocateVector( IntCnts,  j )
!      CALL AllocateVector( IntArray, k )
!
!      Old mesh nodes were copied as is...
!
       IF(.NOT. ASSOCIATED(Mesh % ParallelInfo % Neighbourlist ) ) THEN
         CALL Fatal('UpdateParallelMesh','Original mesh has no NeighbourList!')
       END IF
       
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT. ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
           CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours, 1 )
           NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = ParEnv % MyPe
         ELSE
           CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours, &
               SIZE( Mesh % ParallelInfo % Neighbourlist(i) % Neighbours) )
           NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
               Mesh % ParallelInfo % NeighbourList(i) % Neighbours
         END  IF
       END DO
!
!      Take care of the new mesh internal nodes.
!      Parallel global numbering will take care
!      of the interface nodes:
!      ----------------------------------------
       DO i=Mesh % NumberOfNodes+1, NewMesh % NumberOfNodes
          IF ( .NOT. NewMesh % ParallelInfo % GInterface(i) ) THEN
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
!         IF ( ALL(Mesh % ParallelInfo % GInterface(Edge % NodeIndexes)) .AND. Found ) THEN
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
!            IF ( ALL( Mesh % ParallelInfo % GInterface(Face % NodeIndexes) ) ) THEN
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
       Reorder = [ (i, i=1,NewMesh % NumberOfNodes) ]

       k = NewMesh % Nodes % NumberOfNodes - Mesh % Nodes % NumberOfNodes


       CALL ResetTimer('ParallelGlobalNumbering')
       CALL ParallelGlobalNumbering( NewMesh, Mesh, k, Reorder )
       CALL CheckTimer('ParallelGlobalNumbering',Level=7,Delete=.TRUE.)

       
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
!> Split mesh with quadfaces into one with only triangle faces.
!------------------------------------------------------------------------------
  SUBROUTINE SplitMeshQuads(Mesh, Vlist) 
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Vlist
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    INTEGER :: i, j, k, k2, n, AddCnt, NewElCnt, nBulkElems, nBoundaryElems 
    LOGICAL :: Found, FacesPresent, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Parent
    TYPE(PElementDefs_t), POINTER :: PDefs
    LOGICAL :: Parallel, IsBulkElement, IsOddCut
    LOGICAL :: Is15, Is24, Is35, Is26, Is34, Is16
    INTEGER, ALLOCATABLE :: CutCorner(:), MinCorner(:), BulkElementOffset(:)
    INTEGER :: TypeCnt(0:8), LocalMap(4), LocalMin(3), LocalCut(3)
    INTEGER :: CutChanges, MeshDim, CutComb, iCut
    TYPE(Element_t), POINTER :: Face
    TYPE(Element_t), POINTER :: NewElements(:) => NULL()
    CHARACTER(*), PARAMETER :: Caller="SplitMeshQuads"

!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Mesh ) ) RETURN
    IF(.NOT. ASSOCIATED(VList)) RETURN
    IF(Mesh % MeshDim < 2 ) RETURN
    
    IF( Mesh % MeshDim == 2 ) THEN
      IF(.NOT. ListGetLogical( Vlist,'Split Mesh Quads',Found ) ) RETURN              
    ELSE
      IF(.NOT. ListGetLogical( Vlist,'Split Mesh Prisms',Found ) ) RETURN              
    END IF
    CALL Info( Caller,'Splitting all quadrilaterals into triangles in '//I2S(Mesh % MeshDim)//'D',Level=5)
      
    TypeCnt = 0
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      j = Mesh % Elements(i) % TYPE % ElementCode/100
      TypeCnt(j) = TypeCnt(j) + 1
    END DO
    
    IF(TypeCnt(4) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(4))//' quad elements',Level=8)
    IF(TypeCnt(6) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(6))//' pyramid elements',Level=8)
    IF(TypeCnt(7) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(7))//' prism elements',Level=8)    
    IF(TypeCnt(8) > 0) CALL Info(Caller,'Splitting '//I2S(TypeCnt(8))//' hexahedron elements',Level=8)

    !DO i=0,8
    !  PRINT *,'TypeCount:',i,TypeCnt(i)
    !END DO
    
    IF(TypeCnt(6) + TypeCnt(8) > 0 ) THEN
      CALL Fatal(Caller,'Not implemented yet for pyramids and hexahedrons!')
    END IF
    
    IF(Mesh % MeshDim == 3 .AND. TypeCnt(7) == 0) THEN
      CALL Warn(Caller,'No wedges exist, doing nothing!')
      RETURN
    END IF
    
    CALL ResetTimer(Caller)

    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. Mesh % SingleMesh )
    
    AddCnt = TypeCnt(4) + 2*TypeCnt(7)    
    CALL Info(Caller,'Number of elements added by splitting is '//I2S(AddCnt),Level=6)
    
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z
    MeshDim = Mesh % MeshDim
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    FacesPresent = ASSOCIATED(Mesh % Faces)
        
    NewElCnt = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements + AddCnt
    CALL Info(Caller,'Count of new elements: '//I2S(NewElCnt),Level=7)
    CALL AllocateVector( NewElements, NewElCnt )
    CALL Info(Caller,'New elements allocated.', Level=10 )
    
    ! We do not need to all the complex stuff in 2D.
    IF( MeshDim < 3 ) GOTO 1

    ! First negotiate the split direction of the mesh.
    ! We start from assuming that cut direction is though off local
    ! indexes and the one owning the smallest global index calls the shots.    
    !----------------------------------------------------------------------
    CALL Info(Caller,'Trying to find consistent splitting direction in 3D',Level=12)
    IF(.NOT.FacesPresent) CALL FindMeshFaces3D( Mesh )
    ALLOCATE(CutCorner(Mesh % NumberOfFaces), MinCorner(Mesh % NumberOfFaces))
    
    ! Initialize with smallest index of the face
    DO i=1,Mesh % NumberOfFaces
      Face => Mesh % Faces(i)
      MinCorner(i) = MINVAL(Face % NodeIndexes)
    END DO
    CutCorner = MinCorner

    BLOCK
      INTEGER :: DoMax
      INTEGER, POINTER :: Inds(:)
      REAL(KIND=dp) :: d13,d24
    
      DoMax = 0
      IF(LIstGetLogical(Vlist,'Split Mesh Prisms Min',Found )) DoMax = 1
      IF(LIstGetLogical(Vlist,'Split Mesh Prisms Max',Found )) DoMax = -1

      ! Optionally cut the 3D meshes such that the shorter (longer) diagonal
      ! is used to cut the quad faces. 
      IF(DoMax /= 0) THEN
        DO i=1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i)
          IF(Face % TYPE % ElementCode /= 404) CYCLE
          Inds => Face % NodeIndexes
          
          ! Compute |r1-r3|^2 
          d13 = (x(Inds(1))-x(Inds(3)))**2 + &
              (y(Inds(1))-y(Inds(3)))**2 + (z(Inds(1))-z(Inds(3)))**2 
          
          ! Compute |r2-r4|^2 
          d24 = (x(Inds(2))-x(Inds(4)))**2 + &
              (y(Inds(2))-y(Inds(4)))**2 + (z(Inds(2))-z(Inds(4)))**2 
          
          IF(DoMax * d13 < DoMax * d24) THEN
            CutCorner(i) = MIN(Inds(1),Inds(3))
          ELSE
            CutCorner(i) = MIN(Inds(2),Inds(4))
          END IF
        END DO
      END IF
    END BLOCK
          
    
    
    DO WHILE(.TRUE.)
      CutChanges = 0
      
      DO i=1,Mesh % NumberOfBulkElements         
        Eold => Mesh % Elements(i)        
        SELECT CASE( Eold % TYPE % ElementCode )
          
        CASE(706)
          ! Faces 1 & 2 are triangles
          ! Faces 3, 4 and 5 are married
          ! There are 8 ways to cut the faces but only 6 of those are legal. 
          ! WedgeFaceMap(3,:) = (/ 1,2,5,4 /)
          ! WedgeFaceMap(4,:) = (/ 2,3,6,5 /)
          ! WedgeFaceMap(5,:) = (/ 3,1,4,6 /)
                            
          LocalMin(1:3) = MinCorner(Eold % FaceIndexes(3:5))
          LocalCut(1:3) = CutCorner(Eold % FaceIndexes(3:5))

          ! How are we cutting the three faces? 
          Is24 = ANY( LocalCut(1) == Eold % NodeIndexes([2,4])) 
          Is26 = ANY( LocalCut(2) == Eold % NodeIndexes([2,6])) 
          Is16 = ANY( LocalCut(3) == Eold % NodeIndexes([1,4])) 
          Is15 = .NOT. Is24
          Is35 = .NOT. Is26
          Is34 = .NOT. Is16

          ! Code the cases to numbers [111,222]
          ! Cases 121 and 212 are not allowed!
          CutComb = 111
          IF(.NOT. Is24) CutComb = CutComb+100  ! Is15
          IF(.NOT. Is26) CutComb = CutComb+10   ! Is35
          IF(.NOT. Is16) CutComb = CutComb+1    ! Is34
          
          IF( CutComb == 121 .OR. CutComb == 212 ) THEN
            iCut = 0
            DO k=1,3
              IF(LocalCut(k) == MAXVAL(LocalCut(1:3))) EXIT
            END DO

            IF(k==1) THEN
              IF(CutComb == 121 ) iCut = 1 !-> 221 
              IF(CutComb == 212 ) iCut = 2 !-> 112
            ELSE IF(k==2) THEN
              IF(CutComb == 121 ) iCut = 2 !-> 111
              IF(CutComb == 212 ) iCut = 3 !-> 222
            ELSE IF(k==3) THEN
              IF(CutComb == 121 ) iCut = 3 !-> 122
              IF(CutComb == 212 ) iCut = 1 !-> 211
            ELSE
              CALL Fatal(Caller,'Invalid value for k!')
            END IF

            IF(icut>0) THEN
              CutCorner(Eold % FaceIndexes(2+k)) = Eold % NodeIndexes(iCut)
              CutChanges = CutChanges + 1                            
            ELSE
              CALL Fatal(Caller,'Could not fix invalid cut!')
            END IF
          END IF
        END SELECT
      END DO

      CALL Info(Caller,'Number of switches in cut direction: '//I2S(CutChanges))
      IF(CutChanges == 0) EXIT      
    END DO
    
    ! Jump directly here if we have 2D mesh. 
1   CONTINUE

    ! We need to register the offset coming from split.
    ALLOCATE(BulkElementOffset(Mesh % NumberOfBulkElements))
    BulkElementOffset = 0
    NewElCnt = 0
            
!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      Eold => Mesh % Elements(i)
      IsBulkElement = (i <= Mesh % NumberOfBulkElements ) 
      
      IF(IsBulkElement) THEN
        ! For bulk elements store the offset so we can remap the parents more easily. 
        BulkElementOffset(i) = NewElCnt 
      END IF

      SELECT CASE( Eold % TYPE % ElementCode )

      CASE(101,202,303,504)
        ! Copy elements without quad faces as is. 
        NewElCnt = NewElCnt + 1
        Enew => NewElements(NewElCnt)
        Enew = Eold
        Enew % ElementIndex = NewElCnt         
        CALL AllocateVector( ENew % NodeIndexes, Eold % TYPE % NumberOfNodes )
        Enew % NodeIndexes = Eold % NodeIndexes

        IF(.NOT. IsBulkElement ) THEN
          CALL UpdateParentElements()
        END IF

      CASE(404)
        IF( MeshDim == 3 ) THEN
          Parent => Eold % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED(Parent)) THEN
            Parent => Eold % BoundaryInfo % Right
          END IF          
          Face => Find_Face(Mesh, Eold, Parent ) 
          IsOddCut = ANY( CutCorner(Face % ElementIndex) == Face % NodeIndexes([1,3]) )
        ELSE
          Face => Eold
          IsOddCut = .TRUE.
        END IF
                    
        DO j=1,2         
          IF( IsOddCut ) THEN
            IF(j==1) THEN
              LocalMap(1:3) = [1,2,3]
            ELSE
              LocalMap(1:3) = [3,4,1]
            END IF
          ELSE
            IF(j==1) THEN
              LocalMap(1:3) = [4,1,2]
            ELSE
              LocalMap(1:3) = [2,3,4]
            END IF
          END IF

          NewElCnt = NewElCnt + 1
          Enew => NewElements(NewElCnt)
          Enew = Eold
          Enew % TYPE => GetElementType(303)
          Enew % ElementIndex = NewElCnt         
          CALL AllocateVector( ENew % NodeIndexes, 3 )
          Enew % NodeIndexes(1:3) = Face % NodeIndexes(LocalMap(1:3))
          Enew % EdgeIndexes => NULL()
          Enew % FaceIndexes => NULL()
          
          IF(.NOT. IsBulkElement ) THEN
            CALL UpdateParentElements()
          END IF
        END DO

      CASE(706)

        LocalCut(1:3) = CutCorner(Eold % FaceIndexes(3:5))

        ! How are we cutting the three faces? 
        Is24 = ANY( LocalCut(1) == Eold % NodeIndexes([2,4])) 
        Is26 = ANY( LocalCut(2) == Eold % NodeIndexes([2,6])) 
        Is16 = ANY( LocalCut(3) == Eold % NodeIndexes([1,4])) 
        Is15 = .NOT. Is24
        Is35 = .NOT. Is26
        Is34 = .NOT. Is16

        ! Code the cases to numbers [111,222]
        CutComb = 111
        IF(.NOT. Is24) CutComb = CutComb+100  ! Is15
        IF(.NOT. Is26) CutComb = CutComb+10   ! Is35
        IF(.NOT. Is16) CutComb = CutComb+1    ! Is34
              
        DO j=1,3                   
          SELECT CASE(CutComb)
          CASE(111)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,4,6]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [2,5,4,6]
            ELSE 
              LocalMap(1:4) = [1,2,3,6]
            END IF

          CASE(112)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,4,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [2,5,4,6]
            ELSE 
              LocalMap(1:4) = [2,3,6,4]
            END IF

          CASE(122)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,4,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [2,5,4,3]
            ELSE 
              LocalMap(1:4) = [4,5,3,6]
            END IF

          CASE(221)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,5,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [1,5,4,6]
            ELSE 
              LocalMap(1:4) = [1,3,5,6]
            END IF

          CASE(222)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,5,3]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [1,5,4,3]
            ELSE 
              LocalMap(1:4) = [4,5,3,6]
            END IF

          CASE(211)
            IF(j==1) THEN
              LocalMap(1:4) = [1,2,5,6]
            ELSE IF(j==2) THEN
              LocalMap(1:4) = [1,5,4,6]
            ELSE 
              LocalMap(1:4) = [1,2,3,6]
            END IF

          CASE DEFAULT
            CALL Fatal(Caller,'Unknown case for split: '//I2S(CutComb))

          END SELECT

          NewElCnt = NewElCnt + 1
          Enew => NewElements(NewElCnt)
          Enew = Eold
          Enew % TYPE => GetElementType(504)
          Enew % ElementIndex = NewElCnt         
          CALL AllocateVector( ENew % NodeIndexes, 4 )
          Enew % NodeIndexes = Eold % NodeIndexes(LocalMap(1:4))          
          Enew % EdgeIndexes => NULL()
          Enew % FaceIndexes => NULL()
        END DO
      END SELECT
            
      IF(i==Mesh % NumberOfBulkElements) THEN
        nBulkElems = NewElCnt
      END IF
    END DO

    ! Release old elements and replace them with new elements and element counts
    CALL ReleaseMeshElements( Mesh ) 

    Mesh % Elements => NewElements
    Mesh % NumberOfBulkElements = nBulkElems
    Mesh % NumberOfBoundaryElements = NewElCnt - nBulkElems
            
    ! These are now conservative and could be updated
    ! NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs

    IF( MeshDim == 3 ) THEN
      Mesh % MaxElementNodes = 4
    ELSE
      Mesh % MaxElementNodes = 3
    END IF
    
#if 0
    j = 0
    DO i=1,Mesh % NumberOfBulkElements
      Enew => NewElements(i)        
      IF ( Enew % DGDOFs>0 ) THEN
        Enew % DGDofs == Enew % TYPE % NumberOfNodes
        ALLOCATE(Enew % DGIndexes(Enew % DGDOFs))
        DO k=1,Enew % DGDOFs
          j = j + 1
          Enew % DGIndexes(k)=j
        END DO
      ELSE
        Enew % DGIndexes=>NULL()
      END IF
    END DO
#endif

    DO i=1,NewElCnt 
      IF (i<=Mesh % NumberOfBulkElements) THEN
        Enew => Mesh % Elements(i)
        PDefs => Enew % PDefs
        IF(ASSOCIATED(PDefs)) THEN
          CALL AllocatePDefinitions(Enew)
          Enew % PDefs = PDefs
          
          ! All elements in actual mesh are not edges
          Enew % PDefs % isEdge = .FALSE.
          
          ! If element is of type tetrahedron and is a p element,
          ! do the Ainsworth & Coyle trick
          IF (Enew % TYPE % ElementCode == 504) CALL ConvertToACTetra(Enew)
          CALL GetRefPElementNodes( Enew % TYPE,  Enew % TYPE % NodeU, &
              Enew % TYPE % NodeV, Enew % TYPE % NodeW )
        END IF
      ELSE
        Enew % PDefs=>NULL()
      END IF
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()
      Enew % BubbleIndexes => NULL()
    END DO
    
    CALL CheckTimer(Caller,Delete=.TRUE.)
    
!   Update structures needed for parallel execution:
!   ------------------------------------------------
    IF( Parallel ) THEN
      ! We don't have any updates for parallel mesh.
      ! The nodes stay the same.
    END IF

    ! Release old edges and faces since they don't match the new mesh.
    !-----------------------------------------------------------------
    CALL ReleaseMeshEdgeTables( Mesh )
    CALL ReleaseMeshFaceTables( Mesh )
    !Mesh % Faces => NULL()
    !Mesh % Edges => NULL()
    !Mesh % NumberOfFaces = 0 
    !Mesh % NumberOfEdges = 0 
    
    IF( FacesPresent ) THEN 
      CALL Info(Caller,'Generating faces in the new mesh as they were present in the old!',Level=20)
      CALL FindMeshFaces3D( Mesh )
    END IF
    IF( EdgesPresent ) THEN 
      CALL Info(Caller,'Generating faces in the new mesh as they were present in the old!',Level=20)
      CALL FindMeshEdges( Mesh )
    END IF

    !CALL CheckMeshInfo( Mesh )     
    !CALL writemeshtodisk( Mesh, "koe" )

  CONTAINS

    
    SUBROUTINE UpdateParentElements()

      INTEGER :: j,m,lcnt,l,nCands,ElemCode,BulkOffset
      TYPE(Element_t), POINTER :: CandParent, Parent
      LOGICAL :: hit
      
      ALLOCATE( Enew % BoundaryInfo )
      Enew % BoundaryInfo = Eold % BoundaryInfo 

      DO j=1,2
        IF(j==1) THEN
          Parent => Eold % BoundaryInfo % Left
        ELSE
          Parent => Eold % BoundaryInfo % Right
        END IF
        IF(.NOT. ASSOCIATED( Parent ) ) CYCLE

        ElemCode = Parent % TYPE % ElementCode

        ! Depending on the elementtype the original element has been split into several candidate elements. 
        SELECT CASE(ElemCode)
        CASE(202)
          nCands = 1
        CASE(303)
          nCands = 1
        CASE(404)
          nCands = 2
        CASE(504)
          nCands = 1
        CASE(706)
          nCands = 3
        END SELECT

        hit = .FALSE.

        BulkOffset = BulkElementOffset(Parent % ElementIndex)

        DO k=1,nCands
          CandParent => NewElements(BulkOffset+k)
          m = Enew % Type % NumberOfNodes
          lcnt = 0
          DO l=1,m
            IF(ANY(Enew % NodeIndexes(l) == CandParent % NodeIndexes)) lcnt = lcnt + 1
          END DO
          IF(lcnt == m) THEN
            IF(j==1) THEN
              Enew % BoundaryInfo % Left => CandParent
            ELSE
              Enew % BoundaryInfo % Right => CandParent
            END IF
            hit = .TRUE.
            EXIT
          END IF
        END DO
            
        IF(.NOT. Hit) THEN
          PRINT *,'Not Found:',j,k,lcnt,nCands,ElemCode,Parent % ElementIndex, BulkOffset
          PRINT *,'This:',Eold % NodeIndexes
          PRINT *,'Parent:',Parent % NodeIndexes
          DO k=1,nCands
            CandParent => NewElements(BulkOffset+k)
            PRINT *,'Cands:',CandParent % NodeIndexes
          END DO
          CALL Fatal('UpdateParentElements','Could not find parent for type '//I2S(ElemCode))
        END IF
        
      END DO

    END SUBROUTINE UpdateParentElements
    
     
   END SUBROUTINE SplitMeshQuads
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Sometimes we are lucky and the mesh includes similar elements that are
!> different only by their center point. If we then ensure that their local
!> numbering is the same we may use same finite element basis vectors for them.
!------------------------------------------------------------------------------
  SUBROUTINE SetEqualElementIndeces( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
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

    TYPE(Mesh_t), POINTER :: Mesh
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
    TYPE(Mesh_t), POINTER :: Mesh
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
  SUBROUTINE DetectExtrudedStructure( Mesh, Solver, ExtVar, &
      TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer, &
      MidNodePointer, MidLayerExists, NumberOfLayers, NodeLayer, &
      MaskVar )
    
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
    TYPE(Variable_t), POINTER, OPTIONAL :: MaskVar
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

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
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
    TYPE(Mesh_t), POINTER :: Mesh
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

  SUBROUTINE ElmerMeshToDualGraph(Mesh, DualGraph, UseBoundaryMesh)
    IMPLICIT NONE

    TYPE(Mesh_t) :: Mesh
    TYPE(Graph_t) :: DualGraph
    LOGICAL, OPTIONAL :: UseBoundaryMesh

    TYPE(Element_t), POINTER :: Element, Elements(:)

    ! MESH DATA
    ! Mesh (CRS format)
    INTEGER, ALLOCATABLE :: eptr(:), eind(:)
    INTEGER :: nelem
    ! Vertex to element map (CRS format)
    INTEGER, ALLOCATABLE :: vptr(:), vind(:)
    INTEGER :: nvertex

    ! WORK ARRAYS
    ! Pointers to vertex-element maps of the current element
    INTEGER, ALLOCATABLE :: ptrli(:), ptrti(:)
    ! Neighbour indices
    INTEGER, ALLOCATABLE :: neighind(:)
    ! ARRAY MERGE: map for merge
    INTEGER, ALLOCATABLE :: wrkmap(:)

    TYPE :: IntTuple_t
      INTEGER :: i1, i2
    END type IntTuple_t

    TYPE(IntTuple_t), ALLOCATABLE :: wrkheap(:)

    ! OpenMP thread block leads for work division
    INTEGER, ALLOCATABLE :: thrblk(:)
    ! Work indices
    INTEGER, ALLOCATABLE :: wrkind(:), wrkindresize(:)
    INTEGER :: nwrkind

    ! Variables
    INTEGER :: i, dnnz, eid, nl, nli, nti, nn, nv, nthr, &
            te, thrli, thrti, vli, vti, TID, allocstat
    INTEGER :: mapSizePad, maxNodesPad, neighSizePad
    LOGICAL :: Boundary

    INTEGER, PARAMETER :: HEAPALG_THRESHOLD = 24

    CALL Info('ElmerMeshToDualGraph','Creating a dual graph for the mesh',Level=8)

    Boundary = .FALSE.
    IF (Present(UseBoundaryMesh)) Boundary = UseBoundaryMesh

    ! Pointers to mesh data
    IF (.NOT. Boundary) THEN
       nelem = Mesh % NumberOfBulkElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements
    ELSE
       nelem = Mesh % NumberOfBoundaryElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements(&
            Mesh % NumberOfBulkElements+1:Mesh % NumberOfBulkElements+nelem)
    END IF

    ! Initialize dual mesh size and number of nonzeroes
    DualGraph % n = nelem
    dnnz = 0

    ! Copy mesh to CRS structure
    ALLOCATE(eptr(nelem+1), eind(nelem*Mesh % MaxElementNodes), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate mesh structure!')

    eptr(1)=1 ! Fortran numbering
    DO i=1, nelem
      Element => Elements(i)
      nl = Element % TYPE % NumberOfNodes
      nli = eptr(i) ! Fortran numbering
      nti = nli+nl-1
      eind(nli:nti) = Element % NodeIndexes(1:nl) ! Fortran numbering
      eptr(i+1) = nli+nl
    END DO

    ! Construct vertex to element list (in serial!)
    CALL VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)

    ! Allocate pointers to dual mesh
    ALLOCATE(DualGraph % ptr(nelem+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')

    ! Divide work by number of rows in the vertex graph
    nthr = 1 
    !$ nthr = omp_get_max_threads()

    ! Load balance the actual work done by threads (slow)
    ! CALL ThreadLoadBalanceElementNeighbour(nthr, nelem, eptr, eind, vptr, thrblk)
    CALL ThreadStaticWorkShare(nthr, nelem, thrblk)

    !$OMP PARALLEL SHARED(nelem, nvertex, eptr, eind, &
    !$OMP                 vptr, vind, Mesh, DualGraph, &
    !$OMP                 nthr, thrblk, dnnz) &
    !$OMP PRIVATE(i, eid, nli, nti, nn, nv, vli, vti, te, &
    !$OMP         maxNodesPad, neighSizePad, ptrli, ptrti, &
    !$OMP         wrkheap, wrkmap, neighind, &
    !$OMP         wrkind, nwrkind, wrkindresize, allocstat, &
    !$OMP         mapSizePad, thrli, thrti, TID) NUM_THREADS(nthr) &
    !$OMP DEFAULT(NONE)

    TID = 1
    !$ TID = OMP_GET_THREAD_NUM()+1

    ! Ensure that the vertex to element lists are sorted
    !$OMP DO 
    DO i=1,nvertex
      vli = vptr(i)
      vti = vptr(i+1)-1

      CALL Sort(vti-vli+1, vind(vli:vti))
    END DO
    !$OMP END DO NOWAIT

    ! Allocate work array (local to each thread)
    maxNodesPad = IntegerNBytePad(Mesh % MaxElementNodes, 8)
    neighSizePad = IntegerNBytePad(Mesh % MaxElementNodes*20, 8)

    ! Pointers to vertex maps
    ALLOCATE(neighind(neighSizePad), &
            ptrli(maxNodesPad), ptrti(maxNodesPad), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')
    ! Initialize neighbour indices
    neighind = 0

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      ! With multiple threads, use heap based merge
      ALLOCATE(wrkheap(maxNodesPad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
    ELSE
      ! With a small number of threads, use map -based merge
      mapSizePad = IntegerNBytePad(nelem, 8)
      ALLOCATE(wrkmap(mapSizePad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
      ! Initialize local map
      wrkmap=0
    END IF

    ! Allocate local list for results
    nwrkind = 0
    ALLOCATE(wrkind(nelem/nthr*20), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')

    ! Ensure that all the threads have finished sorting the vertex indices
    !$OMP BARRIER

    ! Get thread indices
    thrli = thrblk(TID)
    thrti = thrblk(TID+1)

    ! For each element
    DO eid=thrli,thrti-1
      nli = eptr(eid)
      nti = eptr(eid+1)-1
      nv = nti-nli+1

      ! Get pointers to vertices related to the nodes of the element
      te = 0
      DO i=nli,nti
        ptrli(i-nli+1)=vptr(eind(i))
        ptrti(i-nli+1)=vptr(eind(i)+1) ! NOTE: This is to make comparison cheaper
        te = te + ptrti(i-nli+1)-ptrli(i-nli+1)
      END DO

      ! Allocate neighind large enough
      IF (SIZE(neighind)<te) THEN
        DEALLOCATE(neighind)
        neighSizePad = IntegerNBytePad(te,8)
        ALLOCATE(neighind(neighSizePad), STAT=allocstat)
        neighind = 0
      END IF

      ! Merge vertex lists (multi-way merge of ordered lists)
      IF (nthr >= HEAPALG_THRESHOLD) THEN
        CALL kWayMergeHeap(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkheap)
      ELSE
        CALL kWayMergeArray(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkmap)
      END IF

      ! Add merged list to final list of vertices
      IF (nn+nwrkind>SIZE(wrkind)) THEN
        ALLOCATE(wrkindresize(MAX(nn+nwrkind,2*SIZE(wrkind))), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                'Unable to allocate local workspace!')
        wrkindresize(1:nwrkind)=wrkind(1:nwrkind)
        DEALLOCATE(wrkind)
        CALL MOVE_ALLOC(wrkindresize, wrkind)
      END IF
      wrkind(nwrkind+1:nwrkind+nn) = neighind(1:nn)
      nwrkind = nwrkind + nn

      ! Store number of row nonzeroes
      DualGraph % ptr(eid)=nn
    END DO

    ! Get the global size of the dual mesh
    !$OMP DO REDUCTION(+:dnnz)
    DO i=1,nthr
      dnnz = nwrkind
    END DO
    !$OMP END DO

    ! Allocate memory for dual mesh indices
    !$OMP SINGLE
    ALLOCATE(DualGraph % ind(dnnz), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')
    ! ptr stores row counts, build crs pointers from them
    CALL ComputeCRSIndexes(nelem, DualGraph % ptr)
    !$OMP END SINGLE

    DualGraph % ind(&
            DualGraph % ptr(thrli):DualGraph % ptr(thrti)-1)=wrkind(1:nwrkind)

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      DEALLOCATE(wrkheap, STAT=allocstat)
    ELSE
      DEALLOCATE(wrkmap, STAT=allocstat)
    END IF
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to deallocate local workspace!')
    DEALLOCATE(neighind, ptrli, ptrti, wrkind)

    !$OMP END PARALLEL

    ! Deallocate the rest of memory
    DEALLOCATE(eind, eptr, vptr, vind, thrblk)

    CALL Info('ElmerMeshToDualGraph','Dual graph created with size '//I2S(dnnz),Level=8)


  CONTAINS

    SUBROUTINE VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nelem, nvertex
      INTEGER :: eptr(:), eind(:)
      INTEGER, ALLOCATABLE :: vptr(:), vind(:)

      INTEGER :: i, j, v, eli, eti, ind, tmpi, tmpip, allocstat

      ! Initialize vertex structure (enough storage for nvertex vertices
      ! having eptr(nelem+1) elements)
      ALLOCATE(vptr(nvertex+1), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
              'Vertex allocation failed!')
      vptr = 0

      ! For each element

      ! Compute number of elements attached to each vertex (size of lists)
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        DO j=eli, eti
          vptr(eind(j))=vptr(eind(j))+1
        END DO
      END DO

      ! Compute in-place cumulative sum (row pointers!)
      CALL ComputeCRSIndexes(nvertex, vptr)

      ! Allocate vertex to element lists
      ALLOCATE(vind(vptr(nvertex+1)), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
              'Vertex allocation failed!')

      ! Construct element lists for each vertex
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        ! For each vertex in element
        DO j=eli, eti
          ! Add connection to vertex eind(j)
          ind = eind(j)
          vind(vptr(ind))=i
          vptr(ind)=vptr(ind)+1
        END DO
      END DO

      ! Correct row pointers
      DO i=nvertex,2,-1
        vptr(i)=vptr(i-1)
      END DO
      vptr(1)=1
    END SUBROUTINE VertexToElementList

    ! k-way merge with an array
    SUBROUTINE kWayMergeArray(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, map)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      INTEGER :: map(:)

      INTEGER :: i, j, k, vindi

      ! Merge nv lists using a map (i.e. an array)
      nn = 1
      DO i=1,nv
        DO j=ptrli(i), ptrti(i)-1
          vindi = vind(j)
          ! Put element to map if it is not already there
          IF (map(vindi)==0 .AND. vindi /= node) THEN
            neighind(nn)=vindi
            ! Increase counter
            map(vindi)=1
            nn=nn+1
          END IF
        END DO
      END DO
      nn=nn-1

      ! Clear map
      DO i=1,nn
        map(neighind(i)) = 0
      END DO
    END SUBROUTINE kWayMergeArray

    ! k-way merge with an actual heap
    SUBROUTINE kWayMergeHeap(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, heap)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      TYPE(IntTuple_t) :: heap(:)

      TYPE(IntTuple_t) :: tmp
      INTEGER :: ii, l, r, mind, ll, tmpval, tmpind

      ! Local variables
      INTEGER :: i, e, nzheap, vindi, lindi, pind

      ! Put elements to heap
      nzheap = 0
      DO i=1,nv
        IF (ptrli(i)<ptrti(i)) THEN
          heap(i) % i1 = vind(ptrli(i))
          heap(i) % i2= i
          ptrli(i) = ptrli(i)+1
          nzheap = nzheap+1
        END IF
      END DO

      ! Build heap
      DO ii=(nzheap/2), 1, -1
        i = ii
        ! CALL BinaryHeapHeapify(heap, nzheap, i)
        DO 
          ! Find index of the minimum element
          IF (2*i<=nzheap) THEN
            IF (heap(2*i) % i1 < heap(i) % i1) THEN
              mind = 2*i
            ELSE
              mind = i
            END IF
            IF (2*i+1<=nzheap) THEN
              IF (heap(2*i+1) % i1 < heap(mind) % i1) mind = 2*i+1
            END IF
          ELSE
            mind = i
          END IF

          IF (mind == i) EXIT

          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO
      END DO

      pind = -1
      nn = 1
      DO e=1,te
        ! Pick the first element from heap
        vindi = heap(1) % i1
        lindi = heap(1) % i2

        ! Remove duplicates
        IF (vindi /= pind .AND. vindi /= node) THEN
          neighind(nn) = vindi
          pind = vindi
          nn = nn+1
        END IF

        ! Add new element from list (if any)
        IF (ptrli(lindi) < ptrti(lindi)) THEN
          heap(1) % i1 = vind(ptrli(lindi))
          heap(1) % i2 = lindi
          ptrli(lindi) = ptrli(lindi)+1
        ELSE
          heap(1) % i1 = heap(nzheap) % i1
          heap(1) % i2 = heap(nzheap) % i2
          nzheap=nzheap-1
        END IF
        ! CALL BinaryHeapHeapify(heap, nzheap, 1)
        i = 1

        DO 
          ! Find the index of the minimum element
          ii = 2*i
          mind = i
          IF (ii+1<=nzheap) THEN
            ! Elements 2*i and 2*i+1 can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
            IF (heap(ii+1) % i1 < heap(mind) % i1) mind = ii+1
          ELSE IF (ii<=nzheap) THEN
            ! Element ii can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
          END IF

          IF (mind == i) EXIT

          ! Bubble down the element
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO

      END DO
      nn=nn-1
    END SUBROUTINE kWayMergeHeap

    SUBROUTINE BinaryHeapHeapify(heap, nelem, sind)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      INTEGER, INTENT(IN) :: sind

      INTEGER :: i, l, r, mind
      TYPE(IntTuple_t) :: tmp

      i = sind
      DO
        l = 2*i
        r = 2*i+1
        ! Find index of the minimum element
        mind = i
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) mind = l
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(mind) % i1) mind = r
        END IF

        IF (mind /= i) THEN
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        ELSE
          EXIT
        END IF
      END DO
    END SUBROUTINE BinaryHeapHeapify

    FUNCTION BinaryHeapIsHeap(heap, nelem) RESULT(heaporder)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      LOGICAL :: heaporder

      INTEGER :: i, l, r

      heaporder = .TRUE.

      DO i=(nelem/2), 1, -1
        l = 2*i
        r = 2*i+1
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'left: ', l, i
            EXIT
          END IF
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'right: ', r, i
            EXIT
          END IF
        END IF
      END DO
    END FUNCTION BinaryHeapIsHeap

  END SUBROUTINE ElmerMeshToDualGraph

  SUBROUTINE Graph_Deallocate(Graph)
    IMPLICIT NONE
    TYPE(Graph_t) :: Graph

    DEALLOCATE(Graph % ptr)
    DEALLOCATE(Graph % ind)
    Graph % n = 0
  END SUBROUTINE Graph_Deallocate

  SUBROUTINE ElmerGraphColour(Graph, Colouring, ConsistentColours)
    IMPLICIT NONE

    TYPE(Graph_t), INTENT(IN) :: Graph
    TYPE(Graphcolour_t) :: Colouring
    LOGICAL, OPTIONAL :: ConsistentColours

    INTEGER, ALLOCATABLE :: uncolored(:)
    INTEGER, ALLOCATABLE :: fc(:), ucptr(:), rc(:), rcnew(:)

    INTEGER :: nc, dualmaxdeg, i, v, w, uci, wci, vli, vti, vcol, wcol, &
            nrc, nunc, nthr, TID, allocstat, gn
    INTEGER, POINTER :: colours(:)
    INTEGER, PARAMETER :: VERTEX_PER_THREAD = 100
    LOGICAL :: consistent

    ! Iterative parallel greedy algorithm (Alg 2.) from 
    ! U. V. Catalyurek, J. Feo, A.H. Gebremedhin, M. Halappanavar, A. Pothen. 
    ! "Graph coloring algorithms for multi-core and massively multithreaded systems".
    ! Parallel computing, 38, 2012, pp. 576--594. 

    ! Initialize number of colours, maximum degree of graph and number of 
    ! uncolored vertices
    nc = 0
    dualmaxdeg = 0
    gn = Graph % n
    nunc = gn

    ! Check if a reproducible colouring is being requested
    consistent = .FALSE.
    IF (PRESENT(ConsistentColours)) consistent = ConsistentColours

    ! Get maximum vertex degree of the given graph
    !$OMP PARALLEL DO SHARED(Graph) &
    !$OMP PRIVATE(v) REDUCTION(max:dualmaxdeg) DEFAULT(NONE)
    DO v=1,Graph % n
      dualmaxdeg = MAX(dualmaxdeg, Graph % ptr(v+1)- Graph % ptr(v))
    END DO
    !$OMP END PARALLEL DO

    nthr = 1
    ! Ensure that each vertex has at most one thread attached to it
    !$ IF (.NOT. consistent) nthr = MIN(omp_get_max_threads(), gn)

    ! Allocate memory for colours of vertices and thread colour pointers
    ALLOCATE(colours(gn), uncolored(gn), ucptr(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate colour maps!')

    !$OMP PARALLEL SHARED(gn, dualmaxdeg, Graph, colours, nunc, &
    !$OMP                 uncolored, ucptr, nthr) &
    !$OMP PRIVATE(uci, vli, vti, v, w, wci, vcol, wcol, fc, nrc, rc, rcnew, &
    !$OMP         allocstat, TID) &
    !$OMP REDUCTION(max:nc) DEFAULT(NONE) NUM_THREADS(nthr)

    TID=1
    !$ TID=OMP_GET_THREAD_NUM()+1

    ! Greedy algorithm colours a given graph with at 
    ! most max_{v\in V} deg(v)+1 colours
    ALLOCATE(fc(dualmaxdeg+1), rc((gn/nthr)+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate local workspace!')
    ! Initialize forbidden colour array (local to thread)
    fc = 0

    ! Initialize colours and uncolored entries
    !$OMP DO 
    DO v=1,gn
      colours(v)=0
      ! U <- V
      uncolored(v)=v
    END DO
    !$OMP END DO

    DO
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1

        ! For each w\in adj(v) do
        DO w=vli, vti
          ! fc[colour[w]]<-v
          !$OMP ATOMIC READ
          wcol = colours(Graph % ind(w))
          IF (wcol /= 0) fc(wcol) = v
        END DO

        ! Find smallest permissible colour for vertex
        ! c <- min\{i>0: fc[i]/=v \}
        DO i=1,dualmaxdeg+1
          IF (fc(i) /= v) THEN
            !$OMP ATOMIC WRITE 
            colours(v) = i
            ! Maintain maximum colour
            nc = MAX(nc, i)
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO

      nrc = 0
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1
        vcol = colours(v)

        ! Make sure that recolour array has enough storage for 
        ! the worst case (all elements need to be added)
        IF (SIZE(rc)<nrc+(vti-vli)+1) THEN
          ALLOCATE(rcnew(MAX(SIZE(rc)*2, nrc+(vti-vli)+1)), STAT=allocstat)
          IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
                  'Unable to allocate local workspace!')
          rcnew(1:nrc)=rc(1:nrc)
          DEALLOCATE(rc)
          CALL MOVE_ALLOC(rcnew, rc)
        END IF

        ! For each w\in adj(v) do
        DO wci=vli,vti
          w = Graph % ind(wci)
          IF (colours(w)==vcol .AND. v>w) THEN
            ! R <- R\bigcup {v} (thread local)
            nrc = nrc + 1
            rc(nrc)=v
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO NOWAIT

      ucptr(TID)=nrc
      !$OMP BARRIER

      !$OMP SINGLE
      CALL ComputeCRSIndexes(nthr, ucptr)
      nunc = ucptr(nthr+1)-1
      !$OMP END SINGLE

      ! U <- R
      uncolored(ucptr(TID):ucptr(TID+1)-1)=rc(1:nrc)
      !$OMP BARRIER

      ! Colour the remaining vertices sequentially if the 
      ! size of the set of uncoloured vertices is small enough
      IF (nunc < nthr*VERTEX_PER_THREAD) THEN
        !$OMP SINGLE
        DO uci=1,nunc
          v = uncolored(uci)
          vli = Graph % ptr(v)
          vti = Graph % ptr(v+1)-1

          ! For each w\in adj(v) do
          DO w=vli, vti
            ! fc[colour[w]]<-v
            wcol = colours(Graph % ind(w))
            IF (wcol /= 0) fc(wcol) = v
          END DO

          ! Find smallest permissible colour for vertex
          ! c <- min\{i>0: fc[i]/=v \}
          DO i=1,dualmaxdeg+1
            IF (fc(i) /= v) THEN
              ! Single thread, no collisions possible 
              colours(v) = i
              ! Maintain maximum colour
              nc = MAX(nc, i)
              EXIT
            END IF
          END DO
        END DO
        !$OMP END SINGLE NOWAIT

        EXIT
      END IF

    END DO

    ! Deallocate thread local storage
    DEALLOCATE(fc, rc)
    !$OMP END PARALLEL

    DEALLOCATE(uncolored, ucptr)

    ! Set up colouring data structure
    Colouring % nc = nc
!   CALL MOVE_ALLOC(colours, Colouring % colours)
    Colouring % colours => colours
  END SUBROUTINE ElmerGraphColour

  SUBROUTINE Colouring_Deallocate(Colours)
    IMPLICIT NONE
    TYPE(GraphColour_t) :: Colours

    DEALLOCATE(Colours % colours)
    Colours % nc = 0
  END SUBROUTINE Colouring_Deallocate

  SUBROUTINE ElmerColouringToGraph(Colours, PackedList)
    IMPLICIT NONE

    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(Graph_t) :: PackedList

    INTEGER, ALLOCATABLE :: cptr(:), cind(:)

    INTEGER :: nc, c, i, n, allocstat

    nc = Colours % nc
    n = size(Colours % colours)
    ALLOCATE(cptr(nc+1), cind(n), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerGatherColourLists','Memory allocation failed.')
    cptr = 0
    ! Count number of elements in each colour
    DO i=1,n
      cptr(Colours % colours(i))=cptr(Colours % colours(i))+1
    END DO

    CALL ComputeCRSIndexes(nc, cptr)

    DO i=1,n
      c=Colours % colours(i)
      cind(cptr(c))=i
      cptr(c)=cptr(c)+1
    END DO

    DO i=nc,2,-1
      cptr(i)=cptr(i-1)
    END DO
    cptr(1)=1

    ! Set up graph data structure
    PackedList % n = nc
    CALL MOVE_ALLOC(cptr, PackedList % ptr)
    CALL MOVE_ALLOC(cind, PackedList % ind)
  END SUBROUTINE ElmerColouringToGraph

  ! Routine constructs colouring for boundary mesh based on colours of main mesh
  SUBROUTINE ElmerBoundaryGraphColour(Mesh, Colours, BoundaryColours)
    IMPLICIT NONE

    TYPE(Mesh_t), INTENT(IN) :: Mesh
    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(GraphColour_t) :: BoundaryColours

    TYPE(Element_t), POINTER :: Element
    INTEGER :: elem, nelem, nbelem, astat, lcolour, rcolour, nbc
    INTEGER, POINTER :: bcolours(:)

    nelem = Mesh % NumberOfBulkElements
    nbelem = Mesh % NumberOfBoundaryElements

    ! Allocate boundary colouring
    ALLOCATE(bcolours(nbelem), STAT=astat)
    IF (astat /= 0) THEN
       CALL Fatal('ElmerBoundaryGraphColour','Unable to allocate boundary colouring')
    END IF
    
    nbc = 0
    ! Loop over boundary mesh
    !$OMP PARALLEL DO &
    !$OMP SHARED(Mesh, nelem, nbelem, Colours, bcolours) &
    !$OMP PRIVATE(Element, lcolour, rcolour) &
    !$OMP REDUCTION(max:nbc) &
    !$OMP DEFAULT(NONE)
    DO elem=1,nbelem       
       Element => Mesh % Elements(nelem+elem)

       ! Try to find colour for boundary element based on left / right parent
       lcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          lcolour = Colours % colours(Element % BoundaryInfo % Left % ElementIndex)
       END IF
       rcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          rcolour = Colours % colours(Element % BoundaryInfo % Right % ElementIndex)
       END IF

       ! Sanity check for debug
       IF (ASSOCIATED(Element % BoundaryInfo % Left) .AND. & 
          ASSOCIATED(Element % BoundaryInfo % Right) .AND. &
            lcolour /= rcolour) THEN
         CALL Warn('ElmerBoundaryGraphColour','Inconsistent colours for boundary element: ' &
               // i2s(elem) // "=>" &
               // i2s(lcolour)// " | "//i2s(rcolour))
         WRITE (*,*) Element % BoundaryInfo % Left % ElementIndex, Element % BoundaryInfo % Right % ElementIndex
       END IF

       bcolours(elem)=MAX(lcolour,rcolour)
       nbc=MAX(nbc,bcolours(elem))
    END DO
    !$OMP END PARALLEL DO

    ! Set up colouring data structure
    BoundaryColours % nc = nbc
!   CALL MOVE_ALLOC(bcolours, BoundaryColours % colours)
    BoundaryColours % colours => bcolours
  END SUBROUTINE ElmerBoundaryGraphColour
  
  ! Given CRS indices, referenced indirectly from graph, 
  ! evenly load balance the work among the nthr threads
  SUBROUTINE ThreadLoadBalanceElementNeighbour(nthr, gn, gptr, gind, &
          rptr, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER :: gptr(:), gind(:), rptr(:)
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, j, k, wrk, gwrk, thrwrk, allocstat

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadLoadBalanceElementNeighbour', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Compute total global work
    gwrk = 0
    DO i=1,gn
      DO j=gptr(i),gptr(i+1)-1
        gwrk = gwrk + (rptr(gind(j)+1)-rptr(gind(j)))
      END DO
    END DO

    ! Amount of work per thread
    thrwrk = CEILING(REAL(gwrk,dp) / nthr)

    ! Find rows for each thread to compute
    blkleads(1)=1
    DO i=1,nthr
      wrk = 0
      ! Acquire enough work for thread i
      DO j=blkleads(i),gn
        DO k=gptr(j),gptr(j+1)-1
          wrk = wrk + (rptr(gind(j)+1)-rptr(gind(j)))
        END DO
        IF (wrk >= thrwrk) EXIT
      END DO

      blkleads(i+1)=j+1
      ! Check if we have run out of rows
      IF (j+1>gn) EXIT
    END DO
    ! Reset number of rows (may be less than or equal to original number)
    nthr = i
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadLoadBalanceElementNeighbour

  SUBROUTINE ThreadStaticWorkShare(nthr, gn, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, rem, thrwrk, allocstat
    INTEGER :: totelem

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadStaticWorkShare', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Assuming even distribution of nodes / element, 
    ! distribute rows for each thread to compute 
    blkleads(1)=1
    thrwrk = gn / nthr
    rem = gn-nthr*thrwrk
    ! totelem = 0
    DO i=1,nthr-1
      IF (i<rem) THEN
        blkleads(i+1)=blkleads(i)+thrwrk+1
      ELSE
        blkleads(i+1)=blkleads(i)+thrwrk
      END IF
    END DO
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadStaticWorkShare

  ! Given row counts, in-place compute CRS indices to data
  SUBROUTINE ComputeCRSIndexes(n, arr)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER :: arr(:)

    INTEGER :: i, indi, indip

    indi = arr(1)
    arr(1)=1
    DO i=1,n-1
      indip=arr(i+1)
      arr(i+1)=arr(i)+indi
      indi=indip
    END DO
    arr(n+1)=arr(n)+indi
  END SUBROUTINE ComputeCRSIndexes

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
    TYPE(Mesh_t), POINTER :: Mesh
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
  FUNCTION SplitMeshLevelset(Mesh,Vlist) RESULT( NewMesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Vlist    
    TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE :: phi(:)
    INTEGER, ALLOCATABLE :: EdgeSplit(:)
    LOGICAL, ALLOCATABLE :: CutNode(:)
    TYPE(Variable_t), POINTER :: Var
    LOGICAL :: SplitReady    
    REAL(KIND=dp), POINTER :: u(:),v(:),w(:),x(:),y(:),z(:)
    REAL(KIND=dp) :: Eps
    INTEGER, POINTER :: NodeIndexes(:), EdgeIndexes(:)    
    INTEGER :: i, j, j2, j3, k, k2, k3, l, l2, l3, m, n, &
        n_old, n_new, n_cut, n_split, n_neg, n_pos
    INTEGER :: NoHits, NewElCnt, BCCnt, prevl, &
        NodeCnt, FaceCnt, Node, ParentId 
    LOGICAL :: Found, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Eptr,Parent 
    INTEGER, POINTER :: Child(:,:)
    REAL(KIND=dp) :: h1,h2,hprod,r,s1,s2 
    REAL(KIND=dp), POINTER :: stime(:)
    INTEGER :: ierr, ParTmp(6), ParSizes(6)
    INTEGER :: BodyOffset, SgnNode, BodyCount, LevelsetBC
    LOGICAL :: PosOffset, BulkParent, Parallel
    CHARACTER(:), ALLOCATABLE :: str       
    CHARACTER(*), PARAMETER :: Caller = 'SplitMeshLevelset'
        
!------------------------------------------------------------------------------
    CALL Info( Caller, 'Splitting finite element mesh at zero levelset!', Level = 5 )

    IF ( .NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Warn(Caller,'Original mesh not associated!')
      RETURN
    END IF
        
    CALL ResetTimer(Caller)
    
    DO i=1,Mesh % NumberOfBulkElements
      n = Mesh % Elements(i) % TYPE % ElementCode 
      IF( n /= 303 .AND. n /= 504 ) THEN
        CALL Fatal(Caller,'Only linear triangles and tets can be split: '//I2S(n))
      END IF
    END DO
    
    Parallel = ( ParEnv % PEs > 1 ) .AND. (.NOT. Mesh % SingleMesh )
        
    CALL Info( Caller, '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( Caller, Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( Caller, Message, Level=6 )

    ! At this stage the coordinates have not been added as variable.
    ! We cannot use the UDF's if these are not available. Also time
    ! is needed by default in some calls. 
    Var => VariableGet( Mesh % Variables,'time')
    IF(.NOT. ASSOCIATED( Var ) ) THEN
      CALL VariableAdd( Mesh % Variables, Mesh, &
          Name='Coordinate 1',DOFs=1,Values=Mesh % Nodes % x )   
      CALL VariableAdd(Mesh % Variables,Mesh, &
          Name='Coordinate 2',DOFs=1,Values=Mesh % Nodes % y )    
      CALL VariableAdd(Mesh % Variables,Mesh, &
          Name='Coordinate 3',DOFs=1,Values=Mesh % Nodes % z )    
      ALLOCATE(stime(1)); stime(1) = 0.0_dp
      CALL VariableAdd( Mesh % Variables, Mesh, &
          Name='Time',DOFs=1, Values=sTime )
      CurrentModel % Variables => Mesh % Variables
    END IF
    
    ! Initialize the levelset function for all nodes
    n_old = Mesh % NumberOfNodes
    ALLOCATE( Phi(n_old) )
    
    str = ListGetString( Vlist,'Levelset Variable', Found)
    IF( Found ) THEN
      Var => VariableGet(Mesh % Variables, str)
      IF(.NOT. ASSOCIATED(Var) ) THEN
        CALL Fatal(Caller,'"Levelset Variable" requested, but not available: '//TRIM(str))
      END IF
      Phi = 1.0_dp
      ! We revert to nodal indexes since it will be easier in the future!
      DO i=1,n_old
        j = Var % Perm(i)
        IF(j>0) Phi(i) = Var % Values(j)
      END DO
    ELSE      
      DO i=1,n_old
        Phi(i) = ListGetRealAtNode(Vlist,'Levelset Function', i, Found)
        IF(.NOT. Found ) THEN
          CALL Fatal(Caller,'"Levelset Function" needed to enrich the mesh!')             
        END IF
      END DO
    END IF
    
    Eps = ListGetCReal( Vlist,'Levelset Epsilon',Found )
    IF(.NOT. Found ) Eps = 1.0e-3
    
    n_pos = COUNT( Phi > 0.0 )
    n_neg = COUNT( Phi < 0.0 ) 
        
    BodyOffset = ListGetInteger( Vlist,'Levelset Body Offset',Found ) 
    PosOffset = ListGetLogical( Vlist,'Levelset Offset Positive',Found ) 
    LevelsetBC = ListGetInteger( Vlist,'Levelset Boundary',Found )
    IF(.NOT. Found) LevelsetBC = CurrentModel % NumberOfBCs
    
    IF( Parallel ) THEN
      n_pos = ParallelReduction(n_pos) 
      n_neg = ParallelReduction(n_neg)
    END IF
    
    CALL Info(Caller,'Positive and negative values: '&
        //I2S(n_pos)//' vs. '//I2S(n_neg),Level=7)    
    
    IF( n_pos == 0 .OR. n_neg == 0 ) THEN
      CALL Warn(Caller,'Nothing to do, no zero levelset available!')
      RETURN
    END IF
       
    ! We need edges in order to do the splitting!
    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT. EdgesPresent) CALL FindMeshEdges( Mesh )
        
    ALLOCATE( EdgeSplit(Mesh % NumberOfEdges), CutNode(n_old) )    
    EdgeSplit = 0
    CutNode = .FALSE.
        
    j = 0
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      h1 = Phi(NodeIndexes(1))
      h2 = Phi(NodeIndexes(2))
      hprod = h1*h2
      IF( hprod < 0.0_dp ) THEN
        r = ABS(h2)/(ABS(h1)+ABS(h2))
        IF( r <= Eps ) THEN
          CutNode(NodeIndexes(2)) = .TRUE.
        ELSE IF(1.0-r < Eps ) THEN
          CutNode(NodeIndexes(1)) = .TRUE.
        ELSE
          j = j+1 
          EdgeSplit(i) = j
        END IF
      ELSE IF( ABS(hprod) < 1.0d-20 ) THEN
        IF(ABS(h1) < 1.0e-20) CutNode(NodeIndexes(1)) = .TRUE. 
        IF(ABS(h2) < 1.0e-20) CutNode(NodeIndexes(2)) = .TRUE.
      END IF
    END DO
    
    n_new = j
    CALL Info(Caller,'Number of additional nodes: '//I2S(n_new),Level=6)

    j = COUNT( CutNode )
    CALL Info(Caller,'Number of cut nodes: '//I2S(j),Level=6)
    
!   Update nodal coordinates:
!   -------------------------
    NodeCnt = n_old + n_new 

!   Create the new mesh
!   -------------------------------
    NewMesh => AllocateMesh()    
    NewMesh % SingleMesh = Mesh % SingleMesh
    NewMesh % Name = Mesh % Name   

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
    x(1:n_old) = u
    y(1:n_old) = v
    z(1:n_old) = w

!   add new nodes where edges are split:
!   ------------------------------------
    DO i=1, Mesh % NumberOfEdges
      NodeIndexes => Mesh % Edges(i) % NodeIndexes
      j = EdgeSplit(i)
      IF( j > 0 ) THEN
        j = j + n_old
        h1 = Phi(NodeIndexes(1))
        h2 = Phi(NodeIndexes(2))
        r = ABS(h2)/(ABS(h1)+ABS(h2))
        x(j) = r*u(NodeIndexes(1)) + (1-r)*u(NodeIndexes(2))
        y(j) = r*v(NodeIndexes(1)) + (1-r)*v(NodeIndexes(2))
        z(j) = r*w(NodeIndexes(1)) + (1-r)*w(NodeIndexes(2))
      END IF
    END DO

    CALL Info(Caller,'Added new nodes on the splitted edges.', Level=10 )  

    
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MinEdgeDOFs = Mesh % MinEdgeDOFs
    NewMesh % MinFaceDOFs = Mesh % MinFaceDOFs
    NewMesh % MaxEdgeDOFs = Mesh % MaxEdgeDOFs
    NewMesh % MaxFaceDOFs = Mesh % MaxFaceDOFs
    NewMesh % MaxElementDOFs = Mesh % MaxElementDOFs
    NewMesh % MeshDim = Mesh % MeshDim

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt

!   Update bulk elements:
!   =====================
!
!   First count maximum number of new elements:
!   -------------------------------------------
    NewElCnt = 0
    BodyCount = 0
    DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Eold => Mesh % Elements(i)
      j = 1

      Found = .FALSE.
      IF( ASSOCIATED( Eold % EdgeIndexes ) ) THEN
        Found = ANY(EdgeSplit(Eold % EdgeIndexes) > 0 )
      ELSE
        CALL Fatal(Caller,'No edges for element: '//I2S(i))
      END IF
      
      IF( Found ) THEN
        SELECT CASE( Eold % TYPE % ElementCode/100 )                
        CASE(2)
          j = 2
        CASE(3)
          j = 3
        CASE(5)
          j = 6
        END SELECT
        ! There will also be additional BC elements on the cut!
        j = j + 1
      END IF
      NewElCnt = NewElCnt + j
    END DO
    
    CALL Info( Caller,'Maximum estimated count of new elements: '//I2S(NewElCnt), Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info(Caller,'New mesh allocated.', Level=20 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 6 )
    Child = 0
    CALL Info(Caller,'Array for bulk elements allocated.', Level=20 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes

!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)
       NodeIndexes => Eold % NodeIndexes       
       n = Eold % TYPE % NumberOfNodes                
       n_split = COUNT( EdgeSplit(Eold % EdgeIndexes) > 0 )

       ! We continue splitting until the element is exhausted
       SplitReady = .FALSE.

       ! Split elements to no more than 6 pieces
       DO m = 1,6
         NewElCnt = NewElCnt + 1
         Child(i,m) = NewElCnt
         Enew => NewMesh % Elements(NewElCnt)

         Enew = Eold
         Enew % TYPE => Eold % TYPE
         Enew % BodyId = Eold % BodyId
         Enew % PartIndex = Eold % PartIndex
         Enew % ElementIndex = NewElCnt
         Enew % NDOFs = Eold % NDOFs
         Enew % EdgeIndexes => NULL()
         Enew % FaceIndexes => NULL()
         Enew % BoundaryInfo => NULL()
         
         CALL AllocateVector( ENew % NodeIndexes, n)
         
         IF( n_split == 0 ) THEN
           Enew % NodeIndexes = NodeIndexes
           DO j=1,n
             IF(.NOT. CutNode(NodeIndexes(j)) ) THEN
               ! This is a representative node that is used to determine the sign of the
               ! new elements in order to decide whether to add offset for body or not. 
               SgnNode = j
               EXIT
             END IF
           END DO
           
           SplitReady = .TRUE.
         ELSE           
           n_cut = COUNT( CutNode(NodeIndexes) )
       
           IF ( Eold % TYPE % ElementCode == 303 ) THEN         
             ! Split triangle to four triangles split on one or two edges
             !-----------------------------------------------------------
             IF( n_split == 2 ) THEN
               DO j=1,3
                 IF( EdgeSplit( Eold % EdgeIndexes(j) ) == 0 ) EXIT
               END DO
               j2 = MODULO(j,3)+1
               j3 = MODULO(j+1,3)+1

               IF( m == 1 ) THEN
                 ! There are two ways to split the triangle.
                 ! Choose the one with shorter diameter.
                 s1 = (x(NodeIndexes(j)) - x(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2 + &
                     (y(NodeIndexes(j)) - y(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2 + &
                     (z(NodeIndexes(j)) - z(n_old + EdgeSplit(Eold % EdgeIndexes(j2))))**2
                 s2 = (x(NodeIndexes(j2)) - x(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2 + &
                     (y(NodeIndexes(j2)) - y(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2 + &
                     (z(NodeIndexes(j2)) - z(n_old + EdgeSplit(Eold % EdgeIndexes(j3))))**2
                 Enew % NodeIndexes(1) = NodeIndexes(j)
                 Enew % NodeIndexes(2) = NodeIndexes(j2)                 
                 IF( s1 < s2 ) THEN
                   Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 ELSE
                   Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                   
                 END IF
                 SgnNode = j
               ELSE IF(m==2) THEN
                 IF( s1 < s2 ) THEN
                   Enew % NodeIndexes(1) = NodeIndexes(j)
                   SgnNode = j
                 ELSE
                   Enew % NodeIndexes(1) = NodeIndexes(j2)                   
                   SgnNode = j2
                 END IF
                 Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 Enew % NodeIndexes(3) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                
               ELSE IF(m==3) THEN
                 Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))
                 Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
                 Enew % NodeIndexes(3) = NodeIndexes(j3)
                 SgnNode = j3
                 SplitReady = .TRUE.
               END IF

             ELSE IF( n_split == 1 ) THEN
               DO j=1,3
                 IF( EdgeSplit( Eold % EdgeIndexes(j) ) > 0 ) EXIT
               END DO
               j2 = MODULO(j,3)+1
               j3 = MODULO(j+1,3)+1

               ! One cut result to splitted elements only if the opposing node is cut through
               IF( .TRUE. .OR. CutNode(NodeIndexes(j3)) ) THEN
                 IF(m==1) THEN
                   Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
                   Enew % NodeIndexes(2) = NodeIndexes(j2)
                   Enew % NodeIndexes(3) = NodeIndexes(j3)
                   IF( CutNode(NodeIndexes(j3)) ) THEN
                     SgnNode = j2
                   ELSE
                     SgnNode = j3
                   END IF
                 ELSE IF(m==2) THEN
                   Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
                   Enew % NodeIndexes(2) = NodeIndexes(j3)
                   Enew % NodeIndexes(3) = NodeIndexes(j)
                   IF( CutNode(NodeIndexes(j3)) ) THEN
                     SgnNode = j
                   ELSE
                     SgnNode = j3
                   END IF 
                   SplitReady = .TRUE.
                 END IF
               ELSE
                 Enew % NodeIndexes = NodeIndexes
                 SgnNode = j3
                 SplitReady = .TRUE.
               END IF
             ELSE
               CALL Fatal(Caller,'Triangle can only deal with 1 and 2 splits!')
             END IF
           ELSE              
             CALL Fatal(Caller,'Element type '//I2S(Eold % TYPE % ElementCode)//&
                 ' not supported by the levelset splitter.')
           END IF
         END IF

         ! Set offset for inside/outside elements of the zero levelset.
         ! The SgnNode is a representative node the sign of which tells whether we are inside
         ! or outside. 
         IF( PosOffset ) THEN
           IF( Phi(NodeIndexes(SgnNode)) > 0.0 )  THEN
             Enew % BodyId = Enew % BodyId + BodyOffset
             BodyCount = BodyCount + 1
           END IF
         ELSE
           IF( Phi(NodeIndexes(SgnNode)) < 0.0 )  THEN
             Enew % BodyId = Enew % BodyId + BodyOffset            
             BodyCount = BodyCount + 1
           END IF
         END IF
         IF( SplitReady ) EXIT
       END DO
     END DO
     
!   Update new mesh element counts:
!   -------------------------------
    NewMesh % NumberOfBulkElements = NewElCnt
    
    CALL Info(Caller,'Number of elements inside: '//I2S(BodyCount),Level=7)
    
   
!   Update boundary elements:
!   ---------------------------------------------------

    BCCnt = 0
    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      IF( i == Mesh % NumberOfBulkElements + 1 ) THEN
        CALL Info(Caller,'Number of boundary elements from bulk cuts: '//I2S(BCCnt))           
        BCCnt = 0
      END IF
     
      Eold => Mesh % Elements(i)
      NodeIndexes => Eold % NodeIndexes             
      BulkParent = ( i <= Mesh % NumberOfBulkElements )
      n_split = COUNT( EdgeSplit(Eold % EdgeIndexes) > 0 )
      n_cut = COUNT( CutNode(NodeIndexes) )

      ! Elements created from bulk cuts require some splits or cuts.
      ! Existing boundary elements remain even without cuts.
      IF( BulkParent ) THEN
        IF( n_split + n_cut <= 1 ) CYCLE
      END IF
  
      SplitReady = .FALSE.
      
      ! Each existing boundary element may be cut to several pieces
      ! For triangles this is just max two!
      DO m=1,10          
        BCCnt = BCCnt + 1
        NewElCnt = NewElCnt + 1
        IF( NewElCnt > SIZE( NewMesh % Elements ) ) THEN
          CALL Fatal(Caller,'Too few elements allocated: '//I2S(NewElCnt))
        END IF
       
        Enew => NewMesh % Elements(NewElCnt)
        
        ALLOCATE(Enew % BoundaryInfo)         
        Enew % PartIndex = Eold % PartIndex
        Enew % ElementIndex = NewElCnt
        
        n = 2
        Enew % TYPE => GetElementType(202)
        CALL AllocateVector( ENew % NodeIndexes, n)
        Enew % NDOFs = n
        Enew % EdgeIndexes => NULL()
        Enew % FaceIndexes => NULL()
                
        IF( BulkParent ) THEN
          ! There are the new boundary elements that come from splitting the mesh
          ! at zero levelset. Give the boundary a new index. 
          Enew % BoundaryInfo % Constraint = LevelsetBC
          
          IF ( Eold % TYPE % ElementCode == 303 ) THEN         
            IF( n_split == 2 ) THEN
              DO j=1,3
                IF( EdgeSplit( Eold % EdgeIndexes(j) ) == 0 ) EXIT
              END DO
              j2 = MODULO(j,3)+1
              j3 = MODULO(j+1,3)+1
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j2))
              Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(j3))                   
            ELSE IF( n_split == 1 .AND. n_cut == 1) THEN
              DO j=1,3
                IF( EdgeSplit( Eold % EdgeIndexes(j) ) > 0 ) EXIT
              END DO
              !j2 = MODULO(j,3)+1
              !j3 = MODULO(j+1,3)+1                        
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(j))
              DO j2=1,3
                IF( CutNode(NodeIndexes(j2)) ) EXIT
              END DO
              Enew % NodeIndexes(2) = NodeIndexes(j2)
            ELSE IF( n_cut == 2) THEN
              DO j=1,3
                IF( .NOT. CutNode(NodeIndexes(j) ) ) EXIT
              END DO
              j2 = MODULO(j,3)+1
              j3 = MODULO(j+1,3)+1                        
              Enew % NodeIndexes(1) = NodeIndexes(j2)
              Enew % NodeIndexes(2) = NodeIndexes(j3)
            ELSE
              CALL Fatal(Caller,'Can only deal with 2 or 1+1 splits!')
            END IF
          ELSE              
            CALL Fatal(Caller,'Element type '//I2S(Eold % TYPE % ElementCode)//&
                ' not supported by the levelset splitting.')
          END IF          
          SplitReady = .TRUE.
          
        ELSE
          ! Each existing boundary element may be cut to several pieces
          Enew % BoundaryInfo = Eold % BoundaryInfo
          
          IF( n_split == 0 ) THEN
            ! If no edge is split the element stays as is
            Enew % NodeIndexes = Eold % NodeIndexes
            SplitReady = .TRUE.
            
          ELSE IF( Eold % TYPE % ElementCode == 202 ) THEN
            IF(m==1) THEN
              Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
              Enew % NodeIndexes(2) = n_old + EdgeSplit(Eold % EdgeIndexes(1))
            ELSE IF(m==2) THEN
              Enew % NodeIndexes(1) = n_old + EdgeSplit(Eold % EdgeIndexes(1))
              Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
              SplitReady = .TRUE.
            END IF
          ELSE
            CALL Fatal(Caller,'Cannot do this element yet!')
          END IF
        END IF

         
        prevl = 0
        DO k=1,2
          ! Pointer to the found left/right bulk element
          Eptr => NULL()

          IF( BulkParent ) THEN
            ! If the boundary results from splitting existing elements then
            ! the parent is the existing bulk elements. 
            Parent => Mesh % Elements(i)
          ELSE
            ! If boundary results from existing boundary elements then the potential
            ! parents are the children of the old parents. 
            IF( k==1 ) THEN
              Parent => Eold % BoundaryInfo % Left
            ELSE            
              Parent => Eold % BoundaryInfo % Right
            END IF
            IF(.NOT. ASSOCIATED(Parent)) CYCLE
          END IF

          ! Find the correct parent among the splitted children of the
          ! initial bulk elements. There may be 1 or several children. 
          DO k2 = 1, 6            
            l = Child( Parent % ElementIndex, k2 )
            IF(l==0) CYCLE
            NoHits = 0
            
            IF( BulkParent ) THEN
              IF( k==2 .AND. l == prevl ) CYCLE
            END IF

            IF(l == 0 .OR. l > SIZE( NewMesh % Elements) ) THEN
            ! This is left for debugging...
#if 1
              PRINT *,'Size Elements:',l, SIZE(NewMesh % Elements)
              PRINT *,'Child:',m,n_split,n_cut,k2,BulkParent
              PRINT *,'Parent index:',i,Parent % ElementIndex
              PRINT *,'Parents children',Child( Parent % ElementIndex, :)
              PRINT *,'Parent indexes:',Parent % NodeIndexes
              PRINT *,'Enew:',Enew % NodeIndexes
              PRINT *,'Eold:',Eold % NodeIndexes
              PRINT *,'Eold edges:',Eold % EdgeIndexes
              PRINT *,'Eold edge node indexes:',Mesh % Edges(Eold % EdgeIndexes(1) ) % NodeIndexes
              DO l2=1,6
                IF(Child( Parent % ElementIndex, l2) == 0 ) EXIT
                PRINT *,'Parent:',l2,NewMesh % &
                    Elements(Child( Parent % ElementIndex,l2)) %  NodeIndexes
              END DO
              PRINT *,'old node indexes:',Mesh % Elements(Parent % ElementIndex) %  NodeIndexes
              PRINT *,'old edge indexes:',Mesh % Elements(Parent % ElementIndex) %  EdgeIndexes
              PRINT *,'old cut indexes:',CutNode(Mesh % Elements(Parent % ElementIndex) %  NodeIndexes)
              PRINT *,'old split indexes:',EdgeSplit(Mesh % Elements(Parent % ElementIndex) %  EdgeIndexes)
#endif
              EXIT
            END IF

            Eptr => NewMesh % Elements(l)

            DO l2 = 1,Enew % Type % NumberOfNodes 
              DO l3 = 1, Eptr % TYPE % NumberOfNodes
                IF( Enew % NodeIndexes(l2) == Eptr % NodeIndexes(l3) ) THEN
                  NoHits = NoHits + 1
                  EXIT
                END IF
              END DO
            END DO
            
            IF( NoHits == n ) EXIT
          END DO
          
          IF( NoHits == n ) THEN
            IF( k==1) THEN
              prevl = l
              Enew % BoundaryInfo % Left => Eptr
            ELSE
              Enew % BoundaryInfo % Right => Eptr
            END IF
          ELSE
            IF(k==1) CALL Warn(Caller,'Could not find even 1 parent!')
          END IF
            
        END DO
       
       ! When we have created all the new boundary elements resulting from splitting
       ! the master element then proceed to next element. 
       IF(SplitReady) EXIT
     END DO
   END DO


!   Update new mesh element counts:
!   -------------------------------
   CALL Info(Caller,'Number of total elements: '//I2S(NewElCnt),Level=7)
    
!   Update new mesh boundary element counts:
!   ----------------------------------------
   NewMesh % NumberOfBoundaryElements = NewElCnt - &
       NewMesh % NumberOfBulkElements
   NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs
   NewMesh % MaxElementNodes = Mesh % MaxElementNodes
   
    
   CALL Info( Caller, '******** New mesh ********', Level=6 )
   WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
   CALL Info( Caller, Message, Level=6 )
   WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
   CALL Info( Caller, Message, Level=6 )
   WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
   CALL Info( Caller, Message, Level=6 )


   ! Information of the new system size, also in parallel
   !----------------------------------------------------------------------
   ParTmp(1) = Mesh % NumberOfNodes
   ParTmp(2) = Mesh % NumberOfBulkElements
   ParTmp(3) = Mesh % NumberOfBoundaryElements
   ParTmp(4) = NewMesh % NumberOfNodes
   ParTmp(5) = NewMesh % NumberOfBulkElements
   ParTmp(6) = NewMesh % NumberOfBoundaryElements
   
   IF( .FALSE. .AND. Parallel ) THEN
     CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
     
     CALL Info(Caller,'Information on parallel mesh sizes',Level=8)
     CALL Info(Caller,'Initial mesh has '//I2S(ParSizes(1))//' nodes',Level=8)
     CALL Info(Caller,'Initial mesh has '//I2S(ParSizes(2))//' bulk elements',Level=8)
     CALL Info(Caller,'Initial mesh has '//I2S(ParSizes(3))//' boundary elements',Level=8)
     CALL Info(Caller,'New mesh has '//I2S(ParSizes(4))//' nodes',Level=5)
     CALL Info(Caller,'New mesh has '//I2S(ParSizes(5))//' bulk elements',Level=5)
     CALL Info(Caller,'New mesh has '//I2S(ParSizes(6))//' boundary elements',Level=5)
   END IF

    
   ! Update structures needed for parallel execution:
   !--------------------------------------------------
   IF( Parallel ) THEN
     CALL UpdateParallelInfo( Mesh, NewMesh )
   END IF

   ! Finalize:
   !-----------
   IF(.NOT.EdgesPresent) THEN
     CALL ReleaseMeshEdgeTables( Mesh )
     CALL ReleaseMeshFaceTables( Mesh )
   ELSE
     CALL FindMeshEdges( NewMesh )
   END IF

   CALL CheckTimer(Caller,Delete=.TRUE.)

   CALL Info(Caller,'Mesh was enriched with zero levelset',Level=8)

 CONTAINS
    
!------------------------------------------------------------------------------
    SUBROUTINE UpdateParallelInfo( Mesh, NewMesh )
!------------------------------------------------------------------------------
      TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
      TYPE(Element_t), POINTER :: Edge
      INTEGER :: i,j1,j2,n,n0,m,istat
      LOGICAL :: Found
!------------------------------------------------------------------------------
!
!      Update mesh interfaces for parallel execution.
!      ==============================================
       n = NewMesh % NumberOfNodes
       ALLOCATE( NewMesh % ParallelInfo % NeighbourList(n), stat=istat )
       IF ( istat /= 0 ) CALL Fatal( Caller, 'Allocate error.' )
       DO i=1,n
         NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       CALL AllocateVector( NewMesh % ParallelInfo % GInterface,n  )       
       NewMesh % ParallelInfo % GInterface = .FALSE.

       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )
       NewMesh % ParallelInfo % GlobalDOFs = 0

       ! Inherit the old parallel data
       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % GInterface(1:n) = Mesh % ParallelInfo % GInterface
       NewMesh % ParallelInfo % GlobalDOFs(1:n) = Mesh % ParallelInfo % GlobalDOFs
       DO i=1,n
         m = SIZE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) 
         ALLOCATE( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours(m) )
         NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
             Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       END DO

       n0 = ParallelReduction(MAXVAL(Mesh % ParallelInfo % GlobalDofs),2)       
       CALL Info(Caller,'Offset for parallel numbering of new nodes: '//I2S(n0))

       ! We need global numbering for the edges that we use for the unique numbering of new nodes
       CALL SParEdgeNumbering(Mesh)
       
       DO i=1,Mesh % NumberOfEdges
         j = EdgeSplit(i)
         IF(j==0) CYCLE
         Edge => Mesh % Edges(j)

         ! Make a unique parallel number for the new nodes introduced at split edges
         NewMesh % ParallelInfo % GlobalDOFs(n+j) = n0 + Edge % GElementIndex         

         j1 = Edge % NodeIndexes(1)
         j2 = Edge % NodeIndexes(2)
         m = CountSameIntegers(Mesh % ParallelInfo % NeighbourList(j1) % Neighbours, &
             Mesh % ParallelInfo % NeighbourList(j2) % Neighbours, &
             NewMesh % ParallelInfo % NeighbourList(n+j) % Neighbours ) 
         NewMesh % ParallelInfo % GInterface(n+j) = (m>1)
       END DO
       
    END SUBROUTINE UpdateParallelInfo
    
  END FUNCTION SplitMeshLevelset
!------------------------------------------------------------------------------


  !> Create interface boundaries consisting of edges defined by the intersection of two higher
  !> dimensional boundaries. This may be useful for 3D meshes where 1D meshes have not been
  !> create in advance.
  !-------------------------------------------------------------------
  SUBROUTINE CreateIntersectionBCs(Model, Mesh)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element, Element2, Enew, Face, Face2, Parent
    INTEGER, POINTER :: NodeIndexes(:), NodeIndexes2(:), EdgeIndexes(:), EdgeIndexes2(:), ParentBCs(:)
    INTEGER :: i,i2,j,j2,k,k2,e,e2,l,n,n2,m,nbc,nbulk,nold,t,t2,istat,newbcs,newcnt,bc_id
    TYPE(Element_t), POINTER :: NewElements(:)
    TYPE(ValueList_t), POINTER :: BC
    INTEGER, ALLOCATABLE :: BoundaryId(:), IntersectionBCs(:,:)
    LOGICAL, ALLOCATABLE :: EdgeDone(:), NodeDone(:)
    LOGICAL :: Found, Hit, EdgesPresent, NeedEdges
    
    ! Count how many of the BCs are intersection BCs that we need to determine
    j = 0
    DO bc_id=1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values
      IF( ListCheckPresent( BC,'Intersection BC' ) .OR. &
          ListCheckPresent( BC,'Intersection Body') ) j = j+1 
    END DO
    NewBCs = j
    IF(NewBCs==0) RETURN

    CALL Info('CreateIntersectionBCs',&
        'Number of intersection BCs to determine: '//I2S(NewBCs),Level=5)

    ! Create a fast look-up table that define the new BC indexes and the parent BCs
    ALLOCATE(IntersectionBCs(j,5))
    IntersectionBCs = 0
    j = 0
    DO bc_id=1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values
      ParentBCs => ListGetIntegerArray( BC,'Intersection BC',Found )
      k = 0 
      IF(.NOT. Found ) THEN
        ! If the intersection is between two bodies mark it separately!
        ParentBCs => ListGetIntegerArray( BC,'Intersection Body',Found )
        k = 1
      END IF
      IF(.NOT. Found) CYCLE
      j = j + 1
      IF(SIZE(ParentBCs) /= 2 ) CALL Fatal('CreateIntersectionBCs','Only available for two parents!')
      IntersectionBCs(j,1) = Model % BCs(bc_id) % Tag
      IntersectionBCs(j,2:3) = ParentBCs(1:2)
      IntersectionBCs(j,4) = k
    END DO

    nbulk = Mesh % NumberOfBulkElements
    nbc = Mesh % NumberOfBoundaryElements
    nold = nbulk + nbc

    ! If we need to find intersection between boundaries create a helper structure. 
    IF( ANY(IntersectionBCs(:,4) == 0) ) THEN
      ALLOCATE( BoundaryId( nbc ) )
      BoundaryId = 0
     
      DO t=1,nbc
        Element => Mesh % Elements(nbulk+t)

        ! Only treat 2D boundary elements
        IF( Mesh % MeshDim == 3 ) THEN
          IF( Element % TYPE % ElementCode < 300 ) CYCLE
        ELSE
          IF( Element % TYPE % ElementCode < 200 ) CYCLE
        END IF

        DO bc_id=1,Model % NumberOfBCs
          IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
        END DO
        IF ( bc_id > Model % NumberOfBCs ) CYCLE

        IF( ANY(IntersectionBCs(:,2)==bc_id .OR. IntersectionBCs(:,3)==bc_id)) THEN
          BoundaryId(t) = bc_id
        END IF
      END DO

      n = COUNT( BoundaryId > 0 )
      CALL Info('CreateIntersectionBCs','Number of candidate intersection parents: '//I2S(n))
    END IF
      
    ! Go the new boundary elements over two times.
    ! On the 1st loop just count the number of new elements.
    ! On the 2nd lopp add the new elements in the element list.
    !-------------------------------------------------------------    
    EdgesPresent = ASSOCIATED( Mesh % Edges )
    NeedEdges = (Mesh % MeshDim == 3 .OR. ANY(IntersectionBCs(:,4)==1) )

    IF(NeedEdges .AND. .NOT. EdgesPresent ) THEN
      CALL Info('CreateInterectionsBCs','Need edges for speedy search!',Level=7)
      CALL FindMeshEdges( Mesh ) 
    END IF

    IF( Mesh % MeshDim == 3 ) THEN
      ALLOCATE( EdgeDone( Mesh % NumberOfEdges ) )       
      CALL CreateIntersection3D(.TRUE.,NewCnt)
      IF(NewCnt==0) THEN
        CALL Info('CreateIntersectionBCs','Could not find any additional interface elements!')        
        GOTO 1
      END IF
      CALL CreateIntersection3D(.FALSE.,NewCnt)
    ELSE
      IF(ANY(IntersectionBCs(:,4) == 0 )) THEN      
        ALLOCATE( NodeDone( Mesh % NumberOfNodes ) )
      END IF
      CALL CreateIntersection2D(.TRUE.,NewCnt)
      IF(NewCnt==0) THEN
        CALL Info('CreateIntersectionBCs','Could not find any additional interface elements!')
        GOTO 1
      END IF
      CALL CreateIntersection2D(.FALSE.,NewCnt)
    END IF
               
    IF( InfoActive(10) ) THEN
      DO i=1,newbcs
        CALL Info('CreateIntersectionBCs','New boundary '//I2S(IntersectionBCs(i,1))//&
            ' with '//I2S(IntersectionBCs(i,5))//' elements')
      END DO
    END IF

1   IF(NeedEdges .AND. .NOT. EdgesPresent ) THEN
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    END IF
      
    CALL Info('CreateIntersectionBCs','All done!',Level=12)

  CONTAINS


    ! Find intersection between 2D boundaries i.e. the result will
    ! be a new 1D boundary, or between 3D bodies and the result will
    ! be a new 2D boundary. 
    !-------------------------------------------------------------
    SUBROUTINE CreateIntersection3D(AllocateOnly,NewCnt)
      LOGICAL :: AllocateOnly
      INTEGER :: NewCnt
      LOGICAL :: BulkMode
            
      EdgeDone = .FALSE.
      NewCnt = 0

      DO l=1,newbcs
        BulkMode = ( IntersectionBCs(l,4) == 1)
        
        IF( BulkMode ) THEN
          DO t=1,Mesh % NumberOfFaces
            Element => Mesh % Faces(t)
            IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

            j = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
              j = Element % BoundaryInfo % Left % BodyId
            END IF
            j2 = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
              j2 = Element % BoundaryInfo % Right % BodyId
            END IF

            IF(j == j2) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            NewCnt = NewCnt + 1

            IF(.NOT. AllocateOnly ) THEN
              Enew => Mesh % Elements(nold+NewCnt)        
              ALLOCATE(Enew % BoundaryInfo)         
              Enew % PartIndex = Element % PartIndex
              Enew % ElementIndex = nold + NewCnt
              Enew % TYPE => GetElementType(Element % Type % ElementCode)

              CALL AllocateVector( ENew % NodeIndexes, n)
              Enew % NodeIndexes = NodeIndexes
              Enew % NDOFs = 1
              Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
              IF(j==0) THEN
                Enew % BoundaryInfo % Left => NULL()
              ELSE
                Enew % BoundaryInfo % Left => Element % BoundaryInfo % Left 
              END IF
              IF(j2==0) THEN
                Enew % BoundaryInfo % Right => NULL()
              ELSE
                Enew % BoundaryInfo % Right => Element % BoundaryInfo % Right 
              END IF

              Enew % EdgeIndexes => NULL()
              Enew % FaceIndexes => NULL()
              Enew % PDefs => NULL()
              Enew % BubbleIndexes => NULL()

              IntersectionBCs(l,5) = IntersectionBCs(l,5) + 1
            END IF
          END DO
        ELSE
          DO t=1,nbc
            j = BoundaryId(t) 
            IF(j==0) CYCLE
            ! Do we have a suitable pair of indexes for the parents
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE

            Element => Mesh % Elements(nbulk+t)
            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            DO t2=t+1,nbc
              j2 = BoundaryId(t2) 
              IF(j2==0) CYCLE
              IF(j==j2) CYCLE
              IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

              Element2 => Mesh % Elements(nbulk+t2)
              NodeIndexes2 => Element2 % NodeIndexes
              n2 = Element2 % TYPE % NumberOfNodes

              ! Do we have any common nodes. Some are required...
              k = 0
              DO i=1,n
                IF( ANY(NodeIndexes(i) == NodeIndexes2(1:n2) ) ) k = k+1
              END DO
              IF(k<2) CYCLE

              EdgeIndexes => Element % EdgeIndexes
              IF(ASSOCIATED(EdgeIndexes)) THEN
                Face => Element
              ELSE
                Face => NULL()
                IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
                  Face => Find_Face(Mesh, Element % BoundaryInfo % Left, Element )                
                END IF
                IF(.NOT. ASSOCIATED(Face) ) THEN
                  IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
                    Face => Find_Face(Mesh, Element % BoundaryInfo % Right, Element )
                  END IF
                END IF
                IF(ASSOCIATED( Face ) ) THEN
                  EdgeIndexes => Face % EdgeIndexes
                ELSE
                  CALL Fatal('CreateIntersectionBCs','EdgeIndexes not associated!')
                END IF
              END IF

              ! This is a probably candidate as we have two 2D elements of proper type
              ! sharing at least two nodes. Just have to find for which edges the intersection
              ! applies. It could sometimes be a false positive also. 
              DO i=1,Face % TYPE % NumberOfEdges
                e = EdgeIndexes(i)          
                IF( EdgeDone(e) ) CYCLE

                EdgeIndexes2 => Element2 % EdgeIndexes
                IF(ASSOCIATED(EdgeIndexes2)) THEN
                  Face2 => Element2
                ELSE
                  Face2 => NULL()
                  IF( ASSOCIATED( Element2 % BoundaryInfo % Left ) ) THEN
                    Face2 => Find_Face(Mesh, Element2 % BoundaryInfo % Left, Element2 )
                  END IF
                  IF(.NOT. ASSOCIATED(Face2) ) THEN
                    IF( ASSOCIATED( Element2 % BoundaryInfo % Right ) ) THEN
                      Face2 => Find_Face(Mesh, Element2 % BoundaryInfo % Right, Element2 )
                    END IF
                  END IF
                  IF(ASSOCIATED( Face2 ) ) THEN
                    EdgeIndexes2 => Face2 % EdgeIndexes
                  ELSE
                    CALL Fatal('CreateIntersectionBCs','EdgeIndexes2 not associated!')
                  END IF
                END IF

                DO i2=1,Face2 % TYPE % NumberOfEdges
                  e2 = EdgeIndexes2(i2)

                  ! Ok, we have a hit. Same edge appearing in the proper parent
                  ! boundary elements. Create the actual boundary element only if
                  ! we have already allocated for it. 
                  IF(e==e2) THEN
                    EdgeDone(e) = .TRUE.
                    NewCnt = NewCnt + 1

                    IF(.NOT. AllocateOnly ) THEN
                      Enew => Mesh % Elements(nold+NewCnt)        
                      ALLOCATE(Enew % BoundaryInfo)         
                      Enew % PartIndex = Element % PartIndex
                      Enew % ElementIndex = nold + NewCnt

                      Enew % TYPE => Mesh % Edges(e) % TYPE

                      m = Enew % TYPE % NumberOfNodes
                      CALL AllocateVector( ENew % NodeIndexes, m)
                      Enew % NodeIndexes = Mesh % Edges(e) % NodeIndexes
                      Enew % NDOFs = m
                      Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
                      Enew % BoundaryInfo % Left => Element
                      Enew % BoundaryInfo % Right => Element2

                      Enew % EdgeIndexes => NULL()
                      Enew % FaceIndexes => NULL()
                      Enew % PDefs => NULL()
                      Enew % BubbleIndexes => NULL()

                      ! Just a simple counter for the new BCs of this type
                      IntersectionBCs(l,4) = IntersectionBCs(l,4) + 1
                    END IF

                    EXIT
                  END IF
                END DO

              END DO
            END DO
          END DO
        END IF
      END DO
        
      ! There is nothing to do since no new elements will be created.
      IF( NewCnt == 0 ) RETURN
              
      IF(AllocateOnly) THEN
        ALLOCATE( NewElements(nold + NewCnt ) )
        CALL Info('CreateIntersectionBCs','Allocated for '//I2S(NewCnt)//' new 1D boundary elements!',Level=6)
        
        NewElements(1:nold) = Mesh % Elements(1:nold)

        DO i=nbulk+1,nold
          Element => Mesh % Elements(i)        
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

          Parent => Element % BoundaryInfo % Left
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
          END IF

          Parent => Element % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
          END IF
        END DO

        DO t=1,Mesh % NumberOfFaces
          Element => Mesh % Faces(t)
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
            Element % BoundaryInfo % Left => &
                NewElements(Element % BoundaryInfo % Left % ElementIndex)
          END IF
          IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
            Element % BoundaryInfo % Right => &
                NewElements(Element % BoundaryInfo % Right % ElementIndex)
          END IF
        END DO
                
        DEALLOCATE(Mesh % Elements)
        Mesh % Elements => NewElements
        Mesh % NumberOfBoundaryElements = nbc + NewCnt
      END IF

    END SUBROUTINE CreateIntersection3D

    
    ! Find intersection between 1D boundaries i.e. the result will
    ! be a new 0D boundary (=node). 
    !-------------------------------------------------------------
    SUBROUTINE CreateIntersection2D(AllocateOnly,NewCnt)
      LOGICAL :: AllocateOnly
      INTEGER :: NewCnt
      LOGICAL :: BulkMode
      
      NewCnt = 0

      DO l=1,newbcs
        BulkMode = ( IntersectionBCs(l,4) == 1)
        
        IF( BulkMode ) THEN
          DO t=1,Mesh % NumberOfEdges
            Element => Mesh % Edges(t)
            IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

            j = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
              j = Element % BoundaryInfo % Left % BodyId
            END IF

            j2 = 0
            IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
              j2 = Element % BoundaryInfo % Right % BodyId
            END IF

            IF(j == j2) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            NewCnt = NewCnt + 1

            IF(.NOT. AllocateOnly ) THEN
              Enew => Mesh % Elements(nold+NewCnt)        
              ALLOCATE(Enew % BoundaryInfo)         
              !Enew % PartIndex = Element % PartIndex
              Enew % ElementIndex = nold + NewCnt
              Enew % TYPE => GetElementType(Element % Type % ElementCode)

              CALL AllocateVector( ENew % NodeIndexes, n)
              Enew % NodeIndexes = NodeIndexes
              Enew % NDOFs = 1
              Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
              IF(j==0) THEN
                Enew % BoundaryInfo % Left => NULL()
              ELSE
                Enew % BoundaryInfo % Left => Element % BoundaryInfo % Left 
              END IF
              IF(j2==0) THEN
                Enew % BoundaryInfo % Right => NULL()
              ELSE
                Enew % BoundaryInfo % Right => Element % BoundaryInfo % Right 
              END IF

              Enew % EdgeIndexes => NULL()
              Enew % FaceIndexes => NULL()
              Enew % PDefs => NULL()
              Enew % BubbleIndexes => NULL()

              IntersectionBCs(l,5) = IntersectionBCs(l,5) + 1
            END IF
          END DO
        ELSE
          DO t=1,nbc
            NodeDone = .FALSE.
            j = BoundaryId(t) 
            IF(j==0) CYCLE
            IF( ALL( IntersectionBCs(l,2:3) /= j)) CYCLE

            Element => Mesh % Elements(nbulk+t)
            NodeIndexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            DO t2=t+1,nbc
              j2 = BoundaryId(t2)
              IF(j2==0) CYCLE

              IF(j==j2) CYCLE
              IF( ALL( IntersectionBCs(l,2:3) /= j2)) CYCLE

              Element2 => Mesh % Elements(nbulk+t2)
              NodeIndexes2 => Element2 % NodeIndexes
              n2 = Element2 % TYPE % NumberOfNodes

              ! Ok, so BC elements t and t2 have suitable indeces. 
              ! Do we have any common nodes. Some are required...
              k = 0
              e = 0
              DO i=1,n
                IF( ANY(NodeIndexes(i) == NodeIndexes2(1:n2) ) ) THEN
                  e = NodeIndexes(i)
                  k = k+1
                END IF
              END DO
              IF(k/=1) CYCLE

              IF(NodeDone(e)) CYCLE
              NodeDone(e) = .TRUE.
              NewCnt = NewCnt + 1
              
              PRINT *,'BC node:',l,t,t2,e,n
              
              IF(.NOT. AllocateOnly ) THEN
                Enew => Mesh % Elements(nold+NewCnt)        
                ALLOCATE(Enew % BoundaryInfo)         
                Enew % PartIndex = Element % PartIndex
                Enew % ElementIndex = nold + NewCnt
                !Enew % TYPE => Element % Type 
                Enew % TYPE => GetElementType(101)

                CALL AllocateVector( ENew % NodeIndexes, 1)
                Enew % NodeIndexes = e
                Enew % NDOFs = 1
                Enew % BoundaryInfo % Constraint = IntersectionBCs(l,1)
                Enew % BoundaryInfo % Left => Element
                Enew % BoundaryInfo % Right => Element2

                Enew % EdgeIndexes => NULL()
                Enew % FaceIndexes => NULL()
                Enew % PDefs => NULL()
                Enew % BubbleIndexes => NULL()
                
                IntersectionBCs(l,5) = IntersectionBCs(l,5) + 1
              END IF
            END DO
          END DO
        END IF
      END DO

      IF( NewCnt == 0 ) RETURN
              
      IF(AllocateOnly) THEN
        ALLOCATE( NewElements(nold + NewCnt ) )
        CALL Info('CreateIntersectionBCs','Allocated for '//I2S(NewCnt)//' new boundary elements in 2D!',Level=6)
        
        NewElements(1:nold) = Mesh % Elements(1:nold)
        
        DO i=nbulk+1,nold
          Element => Mesh % Elements(i)        
          IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE

          Parent => Element % BoundaryInfo % Left
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
          END IF

          Parent => Element % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
          END IF
        END DO

        DO t=1,Mesh % NumberOfEdges
          Element => Mesh % Edges(t)
          IF(ASSOCIATED(Element % BoundaryInfo % Left)) THEN
            Element % BoundaryInfo % Left => &
                NewElements(Element % BoundaryInfo % Left % ElementIndex)
          END IF
          IF(ASSOCIATED(Element % BoundaryInfo % Right)) THEN
            Element % BoundaryInfo % Right => &
                NewElements(Element % BoundaryInfo % Right % ElementIndex)
          END IF
        END DO
        
        DEALLOCATE(Mesh % Elements)
        Mesh % Elements => NewElements
        Mesh % NumberOfBoundaryElements = nbc + NewCnt
      END IF
      
    END SUBROUTINE CreateIntersection2D
    
    
  END SUBROUTINE CreateIntersectionBCs



  !> Sometimes the mesh includes boundaries but it is annoyingly time-consuming to tag
  !> them by hand as their numbers are not known. Then an alternative is to use some simple
  !> rules to detect the existing boundaries and tag them with to the boundary that fulfills
  !> the detection rule.
  !-------------------------------------------------------------------
  SUBROUTINE TagBCsUsingRule(Model, Mesh)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh

    TYPE(Element_t), POINTER :: Element, Parent
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,m,t,t0
    INTEGER :: bc_ind, pSign, rSign, dim, BCsTagged, RuleInd
    TYPE(ValueList_t), POINTER :: BC
    REAL(KIND=dp) :: Coord(3), eps, val, r, rad, phi, phimin, phimax, RuleC
    LOGICAL :: Found, Hit, Parallel, CreateBCs, SplitBC, DoIt
    INTEGER, ALLOCATABLE :: EdgeConstraint(:)
    CHARACTER(:), ALLOCATABLE :: RuleStr
    CHARACTER(*), PARAMETER :: Caller = 'TagBCsUsingRule'
     

    Parallel = ( ParEnv % PEs > 1 )
    IF( Parallel ) THEN
      IF( ListGetLogical( Model % Simulation,'Single Mesh',Found ) ) THEN
        Parallel = .FALSE.
        CALL Info(Caller,'Working on single mesh, so reverting parallel mode to serial!',Level=8)
      END IF
    END IF
    
    ! We may need the rotor radius in defining certain BCs.
    DoIt = ListGetLogical( Model % Simulation,'Rotor Mode',Found) .AND. &
        .NOT. ListCheckPresent( Model % Simulation,'Rotor Radius')
    DoIt = DoIt .OR. ListGetLogical( Model % Simulation,'Determine Rotor Radius',Found)
    IF(DoIt) THEN
      IF( Parallel ) THEN
        CALL Fatal(Caller,'Cannot determine "Rotor Radius" yet in parallel!')
      ELSE        
        Rad = DetermineRotorRadius(Mesh)
        IF(Rad>0) THEN
          CALL ListAddConstReal(Model % Simulation,'Rotor Radius',Rad)
          WRITE(Message,'(A,ES14.6)') '"Rotor Radius" is found to be: ',Rad
          CALL Info(Caller,Message)
        ELSE
          CALL Fatal(Caller,'Could not determine "Rotor Radius", maybe there are not two pieces!?')
        END IF
      END IF
    END IF    
    
    ! Nothing to do with any boundary   
    CreateBCs = ListCheckPresentAnyBC(Model,'Boundary Create')
    IF(.NOT. ( ListCheckPresentAnyBC( Model,'Boundary Levelset') .OR. &
        ListCheckPresentAnyBC(Model,'Boundary Detect') .OR. CreateBCs ) ) RETURN 

    CALL Info(Caller,'Tagging BCs using simple geometric detection rules')

    IF(CreateBCs) THEN
      !IF(.NOT. ASSOCIATED(Mesh % Edges) ) THEN
        CALL FindMeshEdges2D(Mesh)
      !END IF
      ALLOCATE(EdgeConstraint(Mesh % NumberofEdges) )
      EdgeConstraint = 0
    END IF

    SplitBC = .FALSE.
    dim = Mesh % MeshDim
    BCsTagged = 0
    t0 = Mesh % NumberOfBulkElements

    n = 0
    DO t=1, Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t0+t)      
      IF(Element % BoundaryInfo % Constraint == 0) n=n+1
    END DO
    CALL Info(Caller,'Number of unconstrained boundary elements: '//I2S(n))

    IF(.NOT. CreateBCs ) THEN
      m = n
      IF(Parallel) m = ParallelReduction(m)
      IF(m == 0) THEN
        CALL Warn(Caller,'Boundary detection requested but all boundaries already set!')
        RETURN
      END IF
    END IF
    
    n = Mesh % NumberOfNodes 

    DO bc_ind = 1, Model % NumberOfBCs
      BC => Model % BCs(bc_ind) % Values
      
      eps = ListGetCReal(BC,'Boundary Detect Epsilon',Found )
      IF(.NOT. Found) eps = 1.0e-6
    
      pSign = 0
      RuleInd = 0
      !PRINT *,'bc tag:',bc_ind, Model % BCs(bc_ind) % Tag
      
      ! Do we have a simple rule? 
      ! These rules have been designed such that particularly electrical machines
      ! that have rather constant set of BCs can be treated easily. 
      RuleStr = ListGetString(BC,'Boundary Detect',Found )
      IF(.NOT. Found) RuleStr = ListGetString(BC,'Boundary Create',Found )

      ! When we have rules "phimax" and "phimin" they may be augmented with inner and outer
      rSign = 0
      IF( LEN_TRIM(RuleStr) == 12 ) THEN
        IF( RuleStr(8:12) == 'inner' ) THEN
          Rad = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          rSign = -1
        ELSE IF(RuleStr(8:12) == 'outer') THEN
          Rad = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          rSign = 1
        END IF
      END IF

      IF( Found ) THEN
        SELECT CASE( RuleStr )
        CASE('xmin')
          RuleInd = 1
          RuleC = MINVAL( Mesh % Nodes % x )
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('xmax')
          RuleInd = 1
          RuleC = MAXVAL( Mesh % Nodes % x )
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE('ymin')
          RuleInd = 2
          RuleC = MINVAL( Mesh % Nodes % y )
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('ymax')
          RuleInd = 2
          RuleC = MAXVAL( Mesh % Nodes % y )
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE('zmin')
          RuleInd = 3
          RuleC = MINVAL( Mesh % Nodes % z )
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('zmax')
          RuleInd = 3
          RuleC = MAXVAL( Mesh % Nodes % z )
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE('r inner')
          RuleInd = 4 
          RuleC = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          IF(.NOT. Found) CALL Fatal(Caller,'boundary detect "r inner" requires "Rotor Radius" be given!')
          pSign = -1
        CASE('r outer')
          RuleInd = 4
          RuleC = ListGetCReal( CurrentModel % Simulation,'Rotor Radius',Found)
          IF(.NOT. Found) CALL Fatal(Caller,'boundary detect "r outer" requires "Rotor Radius" be given!')
          pSign = 1
        CASE('phimax','phimax inner','phimax outer')
          ! The "inner" and "outer" rules allow to separate rotor and stator
          phimax = -2*PI
          phimin = 2*PI
          m = 0
          DO i=1,n
            ! Skip the nodes at exact origin
            r = SQRT(Mesh % Nodes % x(i)**2 + Mesh % Nodes % y(i)**2)
            IF(r < eps) CYCLE
            IF(rSign == -1 .AND. r > rad-eps ) CYCLE
            IF(rSign == 1 .AND. r < rad+eps ) CYCLE            
            phi = ATAN2(Mesh % Nodes % y(i), Mesh % Nodes % x(i) )
            phimax = MAX(phimax,phi)
            phimin = MIN(phimin,phi)
            m = m+1
          END DO
          IF(Parallel) THEN
            phimin = ParallelReduction(phimin,1)
            phimax = ParallelReduction(phimax,2)
          END IF
          ! There is a discontinuity at PI. Warn about it so we can code some more. 
          IF(phimax - phimin > PI ) CALL Fatal(Caller,'dPhi bigger than PI?')
          RuleC = phimax
          RuleInd = 5           
        CASE('phimin','phimin inner','phimin outer')
          phimax = -2*PI
          phimin = 2*PI
          m  = 0
          DO i=1,n
            r = SQRT(Mesh % Nodes % x(i)**2 + Mesh % Nodes % y(i)**2)
            IF(r < eps) CYCLE
            IF(rSign == -1 .AND. r > rad-eps ) CYCLE
            IF(rSign == 1 .AND. r < rad+eps ) CYCLE            
            phi = ATAN2(Mesh % Nodes % y(i), Mesh % Nodes % x(i) )
            phimax = MAX(phimax,phi)
            phimin = MIN(phimin,phi)
            m = m+1
          END DO
          IF(Parallel) THEN
            phimin = ParallelReduction(phimin,1)
            phimax = ParallelReduction(phimax,2)
          END IF
          IF(phimax - phimin > PI ) CALL Fatal(Caller,'dPhi bigger than PI?')
          RuleC = phimin
          RuleInd = 5          
        CASE('rmin')
          RuleInd = 6
          RuleC = HUGE(RuleC)
          DO i=1,n 
            RuleC = MIN(RuleC, Mesh % Nodes % x(i)**2+Mesh % Nodes % y(i)**2)
          END DO
          RuleC = SQRT(RuleC)
          IF(Parallel) RuleC = ParallelReduction(RuleC,1)
        CASE('rmax')
          RuleInd = 6
          RuleC = 0.0_dp
          DO i=1,n 
            RuleC = MAX(RuleC, Mesh % Nodes % x(i)**2+Mesh % Nodes % y(i)**2)
          END DO
          RuleC = SQRT(RuleC)
          IF(Parallel) RuleC = ParallelReduction(RuleC,2)
        CASE DEFAULT
          CALL Info(Caller,"Available rules: xmin, xmax, ymin, ymax, rmin, rmax, r inner, r outer'&
              //', phimin (inner/outer), phimax (inner/outer)",Level=3)
          CALL Fatal(Caller,'Uknown "Boundary Detect" method: '//TRIM(RuleStr))
        END SELECT
      ELSE          
        IF( .NOT. ListCheckPresent(BC,'Boundary Levelset') ) CYCLE
        ! Should we check the sign of the parent too?
        
        pSign = ListGetInteger(BC,'Boundary Levelset Parent Sign',Found )
        IF(Found) THEN
          IF(ABS(pSign) /= 1 ) THEN
            CALL Fatal(Caller,'"Boundary Levelset Parent Sign" should be either 1 or -1')
          END IF
        END IF        
      END IF

      IF( RuleInd == 3 .AND. Mesh % MeshDim < 3 ) THEN
        CALL Fatal(Caller,'Cannot use z-rules for 2D mesh!')
      END IF
      
      CALL Info(Caller,'Trying to tag elements to boundary: '//I2S(bc_ind),Level=20)
      BCsTagged = 0
      
      CALL TagElements()

      IF(Parallel) BCsTagged = ParallelReduction(BCsTagged)
      
      CALL Info(Caller,'Number of boundary elements "'//TRIM(RuleStr)//'" tagged to '&
          //I2S(bc_ind)//' is: '//I2S(BCsTagged),Level=7)
      IF( BCsTagged == 0 ) THEN
        CALL Fatal(Caller,'Could not find any boundary elements with rule!')
      END IF

    END DO

    IF( CreateBCs ) THEN
      CALL EdgesToBoundaryElements()
      IF(SplitBC) CALL Info(Caller,'Some of the boundaries was an internal one!',Level=10)
    END IF

    CALL Info(Caller,'Done creating additional BCs',Level=10)
    
    
  CONTAINS

    FUNCTION DiffPhi(Coord,Phi0) RESULT( val )
      REAL(KIND=dp) :: Coord(:)
      REAL(KIND=dp) :: Phi0, val
      INTEGER :: i

      IF(SQRT(SUM(Coord(1:2)**2)) < eps) THEN
        val = 0.0_dp
        RETURN
      END IF
      
      val = ATAN2(Coord(2),Coord(1)) - Phi0
      IF( val > PI ) THEN
        val = val - 2*PI
      ELSE IF( val < -PI ) THEN
        val = val + 2*PI
      END IF                 
    END FUNCTION DiffPhi

    
    SUBROUTINE TagElements()      
      INTEGER :: nc, np
      LOGICAL :: Hit
      REAL(KIND=dp) :: val
      INTEGER :: CandElems

      IF(CreateBCs) THEN
        CandElems = Mesh % NumberOfEdges
      ELSE
        CandElems = Mesh % NumberOfBoundaryElements
      END IF
      
      DO t=1, CandElems
        IF(CreateBCs) THEN
          Element => Mesh % Edges(t)
        ELSE
          Element => Mesh % Elements(t0+t)
          IF(Element % BoundaryInfo % Constraint > 0) CYCLE
        END IF
          
        ! Number of corners
        nc = Element % TYPE % ElementCode / 100
        Hit = .TRUE.
        
        DO i=1,nc
          j = Element % NodeIndexes(i)
          Coord(1) = Mesh % Nodes % x(j)
          Coord(2) = Mesh % Nodes % y(j)
          Coord(3) = Mesh % Nodes % z(j)

          SELECT CASE( RuleInd )

          CASE(1,2,3)
            val = Coord(RuleInd) - RuleC
          CASE(4,6)
            val = RuleC - SQRT(SUM(Coord(1:2)**2))                 
          CASE(5)
            val = DiffPhi(Coord,RuleC)
          CASE DEFAULT
            val = ListGetFunVec(BC,'Target Levelset',Coord(1:dim),dim)
          END SELECT
            
          IF(ABS(val) > eps) THEN
            Hit = .FALSE.            
            EXIT
          END IF
        END DO
        IF(.NOT. Hit) CYCLE

        ! We may additionally test inner/outer rule for the radius
        IF(rSign /= 0) THEN
          DO i=1,nc
            j = Element % NodeIndexes(i)
            Coord(1) = Mesh % Nodes % x(j)
            Coord(2) = Mesh % Nodes % y(j)
            Coord(3) = Mesh % Nodes % z(j)
            
            IF( SQRT(SUM(Coord(1:2)**2)) < eps ) CYCLE

            ! This should be negative for rotor, positive for stator
            val = SQRT(SUM(Coord(1:2)**2)) - Rad                
            IF(ABS(val) < eps) CYCLE

            IF( val*rSign < 0 ) THEN
              Hit = .FALSE.            
              EXIT
            END IF
          END DO
          IF(.NOT. Hit) CYCLE
        END IF
          
        
        IF(pSign /= 0) THEN
          IF( ASSOCIATED( Element % BoundaryInfo % Left ) .AND. &
              ASSOCIATED( Element % BoundaryInfo % Right ) ) SplitBC = .TRUE.

          DO k=1,2
            IF(k==1) THEN
              Parent => Element % BoundaryInfo % Left
            ELSE
              Parent => Element % BoundaryInfo % Right
            END IF

            IF(.NOT. ASSOCIATED(Parent)) CYCLE
            
            ! Number of corners, in 3D we must treat tets, pyramids, and wedges.
            np = Parent % TYPE % ElementCode / 100
            IF(np >= 5 .AND. np <= 7) np = np-1
            
            Hit = .TRUE.
            DO i=1,np
              j = Parent % NodeIndexes(i)
              IF(ANY(Element % NodeIndexes(1:nc) == j)) CYCLE
              
              ! Use the 1st non-boundary corner node to test the condition
              Coord(1) = Mesh % Nodes % x(j)
              Coord(2) = Mesh % Nodes % y(j)
              Coord(3) = Mesh % Nodes % z(j)
              
              SELECT CASE( RuleInd )
                
              CASE(4)
                val = RuleC**2 - SUM(Coord(1:2)**2)                 
              CASE DEFAULT
                val = ListGetFunVec(BC,'Target Levelset',Coord(1:dim),dim)
              END SELECT
              
              IF( val*pSign < 0.0_dp ) Hit = .FALSE.
              ! One node is representative
              EXIT
            END DO

            IF(Hit) EXIT
          END DO
          IF(.NOT. Hit) CYCLE
        END IF
        
        IF(CreateBCs) THEN
          EdgeConstraint(t) = Model % BCs(bc_ind) % Tag
        ELSE
          Element % BoundaryInfo % Constraint = Model % BCs(bc_ind) % Tag
        END IF
        BCsTagged = BCsTagged + 1
      END DO
    
    END SUBROUTINE TagElements   

    
    SUBROUTINE EdgesToBoundaryElements()
      
      INTEGER :: nbulk,nbc,nold,nadd,npar,parentcnt(0:2)
      TYPE(Element_t), POINTER :: NewElements(:),Element,Parent,Edge
      
      nbulk = Mesh % NumberOfBulkElements
      nbc = Mesh % NumberOfBoundaryElements
      nold = nbulk + nbc
      nadd = COUNT(EdgeConstraint > 0)

      IF(nadd == 0) THEN
        CALL Info(Caller,'No new boundary elements to add!',Level=6)
        RETURN
      END IF
        
      ALLOCATE( NewElements(nold + nadd ) )
      CALL Info(Caller,'Allocated for '//I2S(nadd)//' new boundary elements!',Level=6)
      
      NewElements(1:nold) = Mesh % Elements(1:nold)
      
      DO i=nbulk+1,nold
        Element => Mesh % Elements(i)        
        IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) CYCLE
        
        Parent => Element % BoundaryInfo % Left
        IF(ASSOCIATED(Parent)) THEN
          NewElements(i) % BoundaryInfo % Left => NewElements(Parent % ElementIndex)
        END IF
        
        Parent => Element % BoundaryInfo % Right
        IF(ASSOCIATED(Parent)) THEN
          NewElements(i) % BoundaryInfo % Right => NewElements(Parent % ElementIndex)
        END IF
      END DO
      
      DEALLOCATE(Mesh % Elements)
      Mesh % Elements => NewElements
      Mesh % NumberOfBoundaryElements = nbc + nadd

      k = nold
      parentCnt = 0

      DO i=1,Mesh % NumberOfEdges
        j = EdgeConstraint(i)
        IF(j==0) CYCLE
        k = k+1

        Element => Mesh % Elements(k)
        Edge => Mesh % Edges(i)

        IF(.NOT. ASSOCIATED(Element)) THEN
          CALL Fatal(Caller,'Element not associated!?')
        END IF
        IF(.NOT. ASSOCIATED(Edge)) THEN
          CALL Fatal(Caller,'Edge not associated!?')
        END IF
          
        Element % TYPE => Edge % TYPE

        IF(.NOT. ASSOCIATED(Element % BoundaryInfo)) THEN
          ALLOCATE(Element % BoundaryInfo)
        END IF

        npar = 0
        IF(ASSOCIATED(Edge % BoundaryInfo) ) THEN
          Parent => Edge % BoundaryInfo % Left        
          IF(ASSOCIATED(Parent)) THEN
            Element % BoundaryInfo % Left => Mesh % Elements(Parent % ElementIndex)
            npar = npar + 1
          END IF
          Parent => Edge % BoundaryInfo % Right
          IF(ASSOCIATED(Parent)) THEN
            Element % BoundaryInfo % Right => Mesh % Elements(Parent % ElementIndex)
            npar = npar + 1
          END IF
        END IF
        Element % BoundaryInfo % Constraint = j

        parentCnt(npar) = parentCnt(npar) + 1
                
        n = Element % TYPE % NumberOfNodes
        ALLOCATE(Element % NodeIndexes(n))
        Element % NodeIndexes = Edge % NodeIndexes
      END DO
      
      DO i=0,2
        j = parentCnt(i)
        IF(j>0) CALL Info(Caller,'New boundary elements with '//I2S(i)//' parents: '//I2S(j),Level=6)
      END DO
      
    END SUBROUTINE EdgesToBoundaryElements
          
  END SUBROUTINE TagBCsUsingRule


  !> Sometimes the mesh includes bodies different from how we would like to resolve
  !> the equations. Then an alternative is to use some simple rule to redefine the
  !> body index.
  !-------------------------------------------------------------------
  SUBROUTINE TagBodiesUsingCondition(Model, Mesh)

    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh

    TYPE(Element_t), POINTER :: Element
    INTEGER :: i,j,k,n,m,t
    INTEGER :: body_id
    REAL(KIND=dp) :: bodyCond(27)
    LOGICAL :: Found, Parallel, Sloppy, Conservative, Switch
    CHARACTER(:), ALLOCATABLE :: str
    CHARACTER(*), PARAMETER :: Caller = 'TagBodiesUsingCondition'

    
    ! Check that there is something to do, else exit
    IF(.NOT. ListCheckPrefix(Model % Simulation,'Body Define Condition') ) RETURN
    
    Parallel = ( ParEnv % PEs > 1 )
    IF( Parallel ) THEN
      IF( ListGetLogical( Model % Simulation,'Single Mesh',Found ) ) THEN
        Parallel = .FALSE.
        CALL Info(Caller,'Working on single mesh, so reverting parallel mode to serial!',Level=8)
      END IF
    END IF
    
    CALL Info(Caller,'Redefining bodies using geometric detection rules')
            
    Sloppy = ListGetLogical( Model % Simulation,'Body Define Sloppy', Found ) 
    Conservative = ListGetLogical( Model % Simulation,'Body Define Conservative', Found ) 

    ! If these are not set we cannot use ListGetReal later on...
    Model % Mesh => Mesh
    Model % Variables => Mesh % Variables
    
    DO i=1,100       
      str = 'Body Define Condition '//I2S(i)
      IF(.NOT. ListCheckPresent(Model % Simulation,str)) EXIT
      
      body_id = ListGetInteger(Model % Simulation,'Body Define Index '//I2S(i),UnfoundFatal=.TRUE.)
      
      k = 0
      DO t=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        Model % CurrentElement => Element
        n = Element % Type % NumberOfNodes

        IF(Element % BodyId == body_id) CYCLE
        
        bodyCond(1:n) = ListGetReal( Model % Simulation, str, n, Element % NodeIndexes )
        m = COUNT(bodyCond(1:n) > 0)
        
        Switch = .FALSE.
        IF( Conservative ) THEN
          Switch = (m==n)
        ELSE IF( Sloppy ) THEN
          Switch = (m>0)
        ELSE
          Switch = (2*m>n)
        END IF

        IF(switch) THEN
          k = k+1
          Element % BodyId = body_id
        END IF
      END DO

      CALL Info(Caller,'Defining body index '//I2S(body_id)//' in '//I2S(k)//' elements')
      
    END DO
      
  END SUBROUTINE TagBodiesUsingCondition
  
!------------------------------------------------------------------------------
END MODULE MeshUtils
!------------------------------------------------------------------------------

!> \}

