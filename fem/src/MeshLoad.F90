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
!>  Mesh loading and preparation: LoadMesh2 and PrepareMesh orchestration.
!>  Extracted from MeshUtils.
!------------------------------------------------------------------------------

MODULE MeshLoad

    USE MeshBasics
    USE MeshSplit
    USE MeshTagging
    IMPLICIT NONE

CONTAINS

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
   Mesh % SingleMesh = (.NOT. Parallel) .AND. (ParEnv % PEs > 1)

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
   CALL LoadMeshStep( 4, BoundariesOnly=BoundariesOnly )

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

 SUBROUTINE PrepareMesh( Model, Mesh, Parallel, Def_Dofs, mySolver )
   TYPE(Model_t) :: Model
   TYPE(Mesh_t), POINTER :: Mesh
   LOGICAL :: Parallel
   INTEGER, OPTIONAL :: Def_Dofs(:,:), mySolver
   TYPE(ValueList_t), POINTER :: Vlist

   LOGICAL :: Found, DoIt
   CHARACTER(*),PARAMETER :: Caller='PrepareMesh'      

   CALL Info(Caller,'Preparing mesh for the future' )

   
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
         UpdateDefDofs(Model % NumberOfSolvers), DG, MultiMesh
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
       MultiMesh = ListGetLogical(Model % Simulation, 'Multiple Meshes', GotIt)
       InDofs = 0
       InDofs(:,1) = 1
       InDofs(:,4) = -1
       DO s=1,Model % NumberOfSolvers
         IF ( MultiMesh ) THEN
           ! Only let a solver contribute its DOF requirements to a mesh it's
           ! actually assigned to -- otherwise every solver's Def_Dofs would
           ! apply to every mesh, over-allocating DOFs on meshes it never
           ! touches. Solver % Mesh is already the resolved mesh pointer, so
           ! this is a direct identity check, no string parsing needed.
           IF ( ASSOCIATED(Model % Solvers(s) % Mesh) ) THEN
             IF ( .NOT. ASSOCIATED(Model % Solvers(s) % Mesh, Mesh) ) CYCLE
           END IF
         END IF
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

       IF (PRESENT(mySolver)) THEN
         DG = Model % Solvers(mySolver) % DG
       ELSE
         DG = .FALSE.
       END IF
       DG = DG .OR. inDofs(el_id,4) == 0
       
       IF ( DG ) THEN
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

       IF (DG) THEN
         NeedEdges = .TRUE.
       ELSE
         ! In the case of a non-DG solver the discontinuous DOFs might be created without
         ! creating edge/face information: 
         NeedEdges = NeedEdges .OR. ANY( inDOFs(el_id,2:3)>0 )
       END IF
!       NeedEdges = NeedEdges .OR. ANY( inDOFs(el_id,2:4)>0 )
       
       
       ! Check if given element is a p element
       IF (FirstOrderElements .AND. inDOFs(el_id,6) > 0) THEN
         CALL AllocatePDefinitions(Element)
         IF (.NOT. DG) NeedEdges = inDOFs(el_id,6) > 1
         
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
       Element % BDOFs = 0
       DO l=1,2        
         IF (l==1) THEN
           Parent => Element % BoundaryInfo % Left
         ELSE
           Parent => Element % BoundaryInfo % Right
         END IF
         IF (.NOT. ASSOCIATED(Parent)) CYCLE

         IF (Parent % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         j = Parent % ElementIndex

         IF (Element % TYPE % DIMENSION == 1) THEN
           IF (j<1 .OR. j>SIZE(EdgeDOFs)) THEN
             IF (ASSOCIATED(Parent % BoundaryInfo)) THEN
               IF (ASSOCIATED(Parent % BoundaryInfo % Left)) THEN
                 j = Parent % BoundaryInfo % Left % ElementIndex
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
           
           IF (.NOT. NeedEdges .AND. InDOFs(el_id,5) > 0) THEN
             IF (.NOT. ASSOCIATED(Element % PDefs)) CALL AllocatePDefinitions(Element)
             Element % PDefs % isEdge = .TRUE.
             Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id,5)))
           ELSE
             Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
           END IF

         END IF
       END DO

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
     IF( Mesh % MaxBDofs > 0 ) THEN
       CALL Info('NonNodalElements','Elementwise (bubble) dofs max: '//I2S(Mesh % MaxBDofs),Level=12)
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

END MODULE MeshLoad
