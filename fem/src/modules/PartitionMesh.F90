!------------------------------------------------------------------------------
!> Partition the finite element mesh using a set of different strategies.
!> This is first a solver, but will perhaps eventually be made an internal routine.
!
!> The idea is that one could use recursive strategies for different pieces of the 
!> finite element mesh. This way the physics can be better taken into account than
!> when using stand alone partitioning tools.
!
!> Currently no parallel operation is supported.
!------------------------------------------------------------------------------
   SUBROUTINE PartitionMeshSolver( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
     TYPE(Model_t), TARGET :: Model    !< All model information (mesh, materials, BCs, etc...)
     REAL(KIND=dp) :: Timestep         !< Timestep size for time dependent simulations
     LOGICAL :: TransientSimulation    !< Steady state or transient simulation
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Params, SectionParams
     INTEGER, ALLOCATABLE :: ParameterInd(:)
     INTEGER, POINTER :: ElementSet(:), ElementPart(:), NodePart(:)
     INTEGER :: NumberOfSets, NumberOfBoundarySets, SetNo, id
     LOGICAL, POINTER :: PartitionCand(:)
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: i,j,k,n, allocstat
     LOGICAL :: Found
     INTEGER, ALLOCATABLE :: EquationPart(:)    

     TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
     CHARACTER(LEN=MAX_NAME_LEN) :: DirectoryName

     !-----------------------------------------------------------------------
     CALL Info('PartitionMesh','Using internal mesh partitioning')


     IF( ParEnv % PEs > 1 ) THEN
       CALL Info('PartitionMes','Currently the mesh can only be partitioned in serial!')
       RETURN
     END IF

     Params => Solver % Values
     Mesh => Solver % Mesh
     n = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
     ALLOCATE( PartitionCand(n), ElementSet(n), ElementPart(n), STAT = allocstat )
     IF( allocstat /= 0 ) THEN
       CALL Fatal('PartitionMesh','Allocation error in partition stuff')
     END IF

     PartitionCand = .FALSE.
     ElementSet = 0
     ElementPart = 0
     
     n = MAX( Model % NumberOfBCs, Model % NumberOfEquations ) 
     ALLOCATE( ParameterInd(n), STAT = allocstat ) 
     IF( allocstat /= 0 ) THEN
       CALL Fatal('PartitionMesh','Allocation error for ParameterInd')
     END IF
     
     ParameterInd = 0
     
     CALL Info('PartitionMesh','Partitioning the boundary elements sets') 

     CALL InitializeBoundaryElementSet(NumberOfBoundarySets)

     DO SetNo = 1, NumberOfBoundarySets
       SectionParams => NULL()
       IF( ParameterInd(SetNo) > 0 ) THEN
         id = ParameterInd(SetNo)
         IF( id <= Model % NumberOfBCs ) THEN
           SectionParams => Model % BCs(id) % Values
         END IF
       END IF
       CALL PartitionMeshPart(SetNo,ElementSet,SectionParams,.TRUE.)
     END DO
     
     IF( NumberOfBoundarySets > 0 ) THEN
       CALL InheritBoundaryPart()
       IF( ListGetLogical( Params,'Partition Mesh Merge Boundaries',Found ) ) THEN
         CALL MergeBoundaryPart()
       END IF
       CALL ExtendBoundaryPart()
     END IF
     
     CALL Info('PartitionMesh','Partition the bulk elements sets')
     CALL InitializeBulkElementSet(NumberOfSets)

     DO SetNo = 1, NumberOfSets
       SectionParams => NULL()
       IF( ParameterInd(SetNo) > 0 ) THEN
         id = ParameterInd(SetNo)
         IF( id <= Model % NumberOfEquations ) THEN
           SectionParams => Model % Equations(id) % Values
         END IF
       END IF
       CALL PartitionMeshPart(SetNo,ElementSet,SectionParams,.FALSE.)
     END DO


     !CALL Info('PartitionMesh','Defining halo mesh')     
     !CALL DefineMeshHalo()
     
     CALL CreateNeighbourList()
     
     DirectoryName = ListGetString( Params,'Output Directory',Found )
     IF( Found ) THEN
       CALL WriteMeshToDiskPartitioned( Model, Mesh,DirectoryName, &
           ElementPart, NeighbourList )
     END IF
     
     IF( ListGetLogical( Params,'Partitioning Exit',Found ) ) STOP
     
     CALL Info('ParitionMesh','Create partitioing variable')
     CALL SetPartitionVariable()

     CALL Info('PartitionMesh','All done for now')


   CONTAINS

     ! Inherit partition from a boundary partition.
     ! In case of conflict the 1st occurrence prevails.
     !-----------------------------------------------------
     SUBROUTINE InheritBoundaryPart()
       
       TYPE(Element_t), POINTER :: Element, Parent
       INTEGER :: t, LeftRight, BoundPart, NoHerited, NoConflict, ElemIndx

       CALL Info('PartitionMesh','Inheriting the boundary patitioning into the bulk mesh') 

       NoHerited = 0
       NoConflict = 0

       DO t=Mesh % NumberOfBulkElements + 1,&
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 
         Element => Solver % Mesh % Elements(t) 
         IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
           DO LeftRight=0,1
             IF( LeftRight == 0 ) THEN
               Parent => Element % BoundaryInfo % Left
             ELSE
               Parent => Element % BoundaryInfo % Right
             END IF

             IF( ASSOCIATED( Parent ) ) THEN
               BoundPart = ElementSet( t ) 
               ElemIndx = Parent % ElementIndex
               IF( ElementSet( ElemIndx ) == 0 ) THEN
                 NoHerited = NoHerited + 1
                 ElementPart( ElemIndx ) = BoundPart
               ELSE IF( ElementSet( ElemIndx ) /= BoundPart ) THEN
                 NoConflict = NoConflict + 1
               END IF
             END IF
           END DO
         END IF
       END DO
       
       CALL Info('PartitionMesh','Number of herited bulk elements: '//TRIM(I2S(NoHerited)))
       CALL Info('PartitionMesh','Number of conflicted bulk elements: '//TRIM(I2S(NoConflict)))
       
     END SUBROUTINE InheritBoundaryPart


     ! Merge partitioned boundaries that share even just one node. 
     ! The algorithm works fine when there is a small number of 
     ! partitions on the boundary and requires minimal additional space.
     !------------------------------------------------------------------
     SUBROUTINE MergeBoundaryPart()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, i, j, MaxPart
       LOGICAL, ALLOCATABLE :: PartFlag(:), PartitionCoupling(:,:)
       INTEGER, ALLOCATABLE :: PartMap(:)
       

       CALL Info('PartitionMesh','Inheriting the boundary patitioning into the bulk mesh') 

       MaxPart = MAXVAL( ElementPart ) 

       ALLOCATE( PartitionCoupling(MaxPart, MaxPart) )
       PartitionCoupling = .FALSE.
       
       ALLOCATE( PartFlag( Mesh % NumberOfNodes ) ) 

       DO i = 1, MaxPart
         PartFlag = .FALSE.
         CALL Info('PartitionMesh','Studying coupling with partition:'//TRIM(I2S(i)),Level=20)
         DO t=1, Mesh % NumberOfBulkElements 
           IF( ElementPart( t ) == i ) THEN 
             Element => Solver % Mesh % Elements(t) 
             PartFlag( Element % NodeIndexes ) = .TRUE.
           END IF
         END DO

         ! Studying to which partitions couple to
         ! Only study cases j>i since the coupling is symmetric
         DO t=1, Mesh % NumberOfBulkElements 
           j = ElementPart( t )
           IF( j > i ) THEN
             Element => Solver % Mesh % Elements(t) 
             IF( ANY( PartFlag( Element % NodeIndexes ) ) ) THEN
               IF( .NOT. PartitionCoupling(i,j) ) THEN
                 CALL Info('PartitionMesh',&
                     'Coupling '//TRIM(I2S(i))//' and '//TRIM(I2S(j)),Level=10)
                 PartitionCoupling(i,j) = .TRUE.
                 PartitionCoupling(j,i) = .TRUE.
               END IF
             END IF
           END IF
         END DO
       END DO

       IF(.NOT. ANY( PartitionCoupling ) ) THEN
         CALL Info('PartitionMesh','Partitions are not coupled')
         RETURN
       END IF

       ! Create mapping of existing partitions to new reduced number of partitions
       ALLOCATE( PartMap( MaxPart ) )
       DO i=1,MaxPart
         PartMap(i) = i
       END DO
       DO i=1,MaxPart
         DO j=i+1,MaxPart
           IF( PartitionCoupling(i,j) ) THEN
             IF( PartMap(i) /= PartMap(j) ) THEN
               CALL Info('PartitionMesh','Mapping partition '&
                   //TRIM(I2S(j))//' to be '//TRIM(I2S(PartMap(i))),Level=8)
               PartMap(j) = PartMap(i)
             END IF
           END IF
         END DO
       END DO
       j = MAXVAL( PartMap ) 
       CALL Info('PartitionMesh','Number of mapped partitions: '//TRIM(I2S(j)))

       ! The coupling is studied via bulk elements as they are all that matters. 
       ! In the end we also remap the boundary elements for consistancy. 
       DO t=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         i = ElementPart( t ) 
         IF( i > 0 ) THEN
           ElementPart(t) = PartMap(i)
         END IF
       END DO
       CALL Info('PartitionMesh','Connected boundaries merged')

       
     END SUBROUTINE MergeBoundaryPart




     ! Extend partition from an existing bulk partitioning. 
     ! In case of conflict the dominating partitioning prevails.
     ! The routine is written with just a small number of existing
     ! boundary partitions in mind and uses minimal memory. 
     !------------------------------------------------------------

     SUBROUTINE ExtendBoundaryPart()
       
       TYPE(Element_t), POINTER :: Element
       INTEGER :: t, ExtendLayers, NoExtend, ElemIndx, NoHits, TestPart, NumberOfParts
       LOGICAL, ALLOCATABLE :: ActiveNode(:)
       INTEGER, ALLOCATABLE :: RefHits(:)
       INTEGER :: allocstat

       NoExtend = 0

       ExtendLayers = ListGetInteger( Params,'Partition Mesh Extend Layers', Found ) 
       IF( ExtendLayers <= 0 ) RETURN
       
       CALL Info('PartitionMesh','Extending boundary meshes by layers: '//TRIM(I2S(ExtendLayers)))

       ALLOCATE( ActiveNode( Mesh % NumberOfNodes ), STAT = allocstat ) 
       IF( allocstat /= 0 ) THEN
         CALL Fatal('PartitionMesh','Allocation error for ActiveNode')
       END IF


       NumberOfParts = MAXVAL( ElementPart ) 
       IF( NumberOfParts > 1 ) THEN
         CALL Info('PartitionMesh','Extending boundary to dominating owner among: '//TRIM(I2S(NumberOfParts)))
         ALLOCATE( RefHits( Mesh % NumberOfBulkElements ), STAT = allocstat )
         IF( allocstat /= 0 ) THEN
           CALL Fatal('PartitionMesh','Allocation error for RefHits')
         END IF
       END IF


       DO i=1,ExtendLayers

         ! If testing for several partitions then nullify the reference
         IF( NumberOfParts > 1 ) THEN
           RefHits = 0
         END IF

         DO TestPart = 1, NumberOfParts         

           ! Set the active nodes for the partition under testing
           ActiveNode = .FALSE.
           DO t=1, Mesh % NumberOfBulkElements 
             IF( ElementPart( t ) == TestPart ) THEN
               Element => Solver % Mesh % Elements(t) 
               ActiveNode( Element % NodeIndexes ) = .TRUE.
             END IF
           END DO
           
           ! Count the number of hits for this partition
           ! If larger than the maximum so far set the partition
           ! For just one existing partitioning no checks need to be done. 
           !--------------------------------------------------------------
           DO t=1, Mesh % NumberOfBulkElements 
             IF( ElementPart( t ) /= 0 ) CYCLE

             Element => Solver % Mesh % Elements(t) 
             NoHits = COUNT( ActiveNode( Element % NodeIndexes ) )
             IF( NoHits == 0 ) CYCLE

             IF( NumberOfParts > 1 ) THEN
               IF( NoHits <= RefHits( t ) ) CYCLE
               RefHits( t ) = NoHits
             END IF
             
             ElementPart( t ) = TestPart 
             NoExtend = NoExtend + 1
           END DO
         END DO
       END DO

       CALL Info('PartitionMesh','Number of extended bulk elements: '//TRIM(I2S(NoExtend)))
       
       DEALLOCATE( ActiveNode ) 
       IF( NumberOfParts > 1 ) DEALLOCATE( RefHits ) 

     END SUBROUTINE ExtendBoundaryPart


     ! Create a variable for the output of the partitioning.
     ! This is an elemental field, not a nodal one. 
     !------------------------------------------------------
     SUBROUTINE SetPartitionVariable()
       
       TYPE(Variable_t), POINTER :: Var
       INTEGER :: t
       INTEGER, POINTER :: DGIndexes(:)
       TYPE(Element_t), POINTER :: Element
       INTEGER :: n

       Var => Solver % Variable 
       IF( .NOT. ASSOCIATED( Var ) ) RETURN
       
       CALL Info('PartitionMesh','Setting discontinuous field > Partition < ')
       
       Var % Values = 0.0_dp
       
       n = Mesh % NumberOfBulkElements 
       DO t=1,n
         Element => Solver % Mesh % Elements(t) 
         IF( ASSOCIATED( Element % DGIndexes ) ) THEN
           Var % Values( Element % DGIndexes ) = 1.0_dp * ElementPart(t)
         ELSE 
           Var % Values( Element % NodeIndexes ) = 1.0_dp * ElementPart(t)
         END IF
       END DO

     END SUBROUTINE SetPartitionVariable
      

     ! Initialize sets of boundary elements to be partitioned with various strategies
     !--------------------------------------------------------------------------
     SUBROUTINE InitializeBoundaryElementSet( NumberOfParts )
       
       INTEGER :: NumberOfParts
       
       INTEGER :: i,j,k,bc_id
       TYPE(ValueList_t), POINTER :: ValueList
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: Found, SeparateBoundarySets        
       INTEGER, ALLOCATABLE :: BCPart(:)
       INTEGER :: allocstat

       
       SeparateBoundarySets = ListGetLogical( Params, &
           'Partitioning Separate Boundary Set', Found)
       
       ALLOCATE( BCPart( Model % NUmberOfBCs ), STAT = allocstat )
       IF( allocstat /= 0 ) THEN
         CALL Fatal('PartitionMesh','Allocation error for BCPart')
       END IF
       BCPart = 0

       NumberOfParts = 0
       DO bc_id = 1, Model % NumberOfBCs 
         ValueList => Model % BCs(bc_id) % Values
         k = ListGetInteger( ValueList,'Partition Set',Found)
         IF( .NOT. Found ) CYCLE
         BCPart(bc_id) = k 
         ParameterInd(k) = bc_id
       END DO
       

       j = 0
       DO WHILE( ANY( BCPart == j+1 ) )
         j = j + 1
       END DO
       IF( j == 0 ) RETURN


       DO bc_id = 1, Model % NumberOfBCs 
         ValueList => Model % BCs(bc_id) % Values
         
         IF( BCPart(bc_id) > 0 ) CYCLE
         
         IF( ListGetLogical( ValueList, 'Discontinuous Boundary', Found) ) THEN
           BCPart(bc_id) = j
           ParameterInd(j) = bc_id
         END IF
         
         k = ListGetInteger( ValueList, 'Periodic Boundary',Found) 
         IF( k > 0 ) THEN
           BCPart(bc_id) = j
           BCPart(k) = j
           ParameterInd(j) = bc_id
         END IF
         
         k = ListGetInteger( ValueList, 'Mortar Boundary',Found) 
         IF( k > 0 ) THEN
           BCPart(bc_id) = j
           BCPart(k) = j
           ParameterInd(j) = bc_id
         END IF
         
         IF( SeparateBoundarySets .AND. BCPart(bc_id) > 0 ) THEN
           DO WHILE( ANY( BCPart == j ) .OR. ANY( EquationPart == j ) ) 
             j = j + 1
           END DO
         END IF
       END DO
      
       NumberOfParts = MAXVAL( BCPart ) 
       
       IF( NumberOfParts == 0 ) RETURN
       

       DO i = Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         Element => Mesh % Elements(i)
         
         DO bc_id=1,Model % NumberOfBCs
           IF ( Element % BoundaryInfo % Constraint == &
               Model % BCs(bc_id) % Tag ) THEN
             IF( BCPart( bc_id ) > 0 ) THEN
               ElementSet( i ) = BCPart( bc_id )
             END IF
             EXIT
           END IF
         END DO
       END DO

       CALL Info('PartitionMesh','Number of sets for boundary partitioning: '//TRIM(I2S(NumberOfParts)))
       CALL Info('PartitionMesh','Number of boundary elements set: '//TRIM(I2S(COUNT(ElementSet > 0))))
      
     END SUBROUTINE InitializeBoundaryElementSet


     
     ! Initialize sets of bulk elements to be partitioned with various strategies
     ! By default there is only one set but certain BCs are treated differently.
     !--------------------------------------------------------------------------
     SUBROUTINE InitializeBulkElementSet( NumberOfParts )
       
       INTEGER :: NumberOfParts
       
       INTEGER :: i,j,k,eq_id, bc_id
       TYPE(ValueList_t), POINTER :: ValueList
       TYPE(Element_t), POINTER :: Element
       LOGICAL :: Found, SeparateBoundarySets, FoundAny 
       INTEGER :: allocstat
     
       ElementSet = 0 
       
       ALLOCATE( EquationPart( Model % NumberOfEquations ), STAT = allocstat ) 
       IF( allocstat /= 0 ) THEN
         CALL Fatal('PartitionMesh','Allocation error for EquationPart')
       END IF
 
       EquationPart = 0
       FoundAny = .FALSE.
       DO eq_id = 1, Model % NumberOfEquations 
         ValueList => Model % Equations(eq_id) % Values
         k = ListGetInteger( ValueList,'Partition Set',Found) 
         IF( k > 0 ) THEN
           EquationPart(eq_id) = k 
           FoundAny = .TRUE.
         END IF
       END DO
       
       NumberOfParts = MAXVAL( EquationPart ) 
      
       Found = .FALSE.
       DO i = 1, Mesh % NumberOfBulkElements 
         IF( ElementPart(i) > 0 ) CYCLE

         Element => Mesh % Elements(i)
         j = 0
         IF( NumberOfParts > 0 ) THEN
           eq_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values, &
               'Equation', Found )
           IF( eq_id > 0 ) j = EquationPart( eq_id ) 
         END IF
         
         IF( j == 0 ) THEN
           Found = .TRUE.
           ElementSet(i) = NumberOfParts + 1
         ELSE
           ElementSet( i ) = j
         END IF
       END DO
       IF( Found ) NumberOfParts = NumberOfParts + 1
       
       CALL Info('PartitionMesh','Number of sets for bulk partitioning: '//TRIM(I2S(NumberOfParts)))
       
     END SUBROUTINE InitializeBulkElementSet


     
     ! Partition the nodes that have the correct ElementSet using various strategies.
     !------------------------------------------------------------------------
     SUBROUTINE PartitionMeshPart(SetNo, ElementSet, LocalParams, IsBoundary )
       
      INTEGER :: SetNo
      INTEGER, POINTER :: ElementSet(:)
      LOGICAL :: IsBoundary 
      TYPE(ValueList_t), POINTER :: LocalParams 

      
      LOGICAL :: BoundaryPart
      CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform, SetMethod
      LOGICAL :: GotCoordTransform, SetNodes
      INTEGER :: SumPartitions, NoPartitions, NoCand
      LOGICAL :: Found
      INTEGER :: NoCandElements
      REAL(KIND=dp) :: BoundaryFraction
      INTEGER :: PartOffset

      
      PartitionCand = ( ElementSet == SetNo )
      n = Solver % Mesh % NumberOfBulkElements
      
      NoCandElements = COUNT( PartitionCand ) 


      CALL Info('PartitionMesh','Doing element set: '//TRIM(I2S(SetNo)))
      CALL Info('PartitionMesh','Number of elements in set: '//TRIM(I2S(NoCandElements)))

      IF( NoCandElements == 0 ) RETURN

      
      BoundaryFraction = ListGetCReal( Params,&
          'Boundary Partitioning Maximum Fraction',Found)
      IF( IsBoundary .AND. NoCandElements <= &
          BoundaryFraction * Solver % Mesh % NumberOfBulkElements ) THEN
        WRITE(Message,'(A,ES12.3)') 'Number of boundary elements below critical limit: ',BoundaryFraction
        CALL Info('PartitionMesh',Message )
        WHERE( PartitionCand ) ElementPart = SetNo
        RETURN
      END IF
         
      BoundaryPart = .FALSE.

      Mesh => Model % Mesh

      Found = .FALSE.
      IF( ASSOCIATED( LocalParams) ) THEN
        SetMethod = ListGetString( LocalParams,'Partitioning Method',Found)
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          SetMethod = ListGetString( Params,'Boundary Partitioning Method',Found)
        END IF
        IF(.NOT. Found ) SetMethod = ListGetString( Params,'Partitioning Method',Found)
      END IF
      IF( Found ) THEN
        CALL Info('PartitionMesh','Using partition method: '//TRIM(SetMethod))
      ELSE
        CALL Fatal('PartitionMesh','Could not define > Partitioning Method < ')
      END IF

      Found = .FALSE.
      IF( ASSOCIATED( LocalParams) ) THEN
        CoordTransform = ListGetString( LocalParams,&
            'Partitioning Coordinate Transformation',Found)
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          CoordTransform = ListGetString( Params,&
              'Boundary Partitioning Coordinate Transformation',Found)
        ELSE
          CoordTransform = ListGetString( Params,&
              'Partitioning Coordinate Transformation',Found)
        END IF
      END IF
      GotCoordTransform = Found

      IF( ASSOCIATED( LocalParams) ) THEN
        SetNodes = ListGetLogical( LocalParams,'Partition Nodes',Found )
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          SetNodes = ListGetLogical( Params,'Boundary Partition Nodes',Found)
        ELSE
          SetNodes = ListGetLogical( Params,'Partition Nodes',Found)
        END IF
      END IF

      IF( ASSOCIATED( LocalParams) ) THEN
        NoPartitions = ListGetInteger( LocalParams,'Number of Partitions',Found )
      END IF
      IF(.NOT. Found ) THEN
        IF( IsBoundary ) THEN
          NoPartitions = ListGetInteger( Params,'Boundary Number of Partitions',Found)
        ELSE
          NoPartitions = ListGetInteger( Params,'Number Of Partitions',Found)
        END IF
      END IF


      ! There may be various coordinate transformation (e.g. to cylindrical coordinates)
      ! that allow for different partitions when using the geometries partitioning routines. 
      IF( GotCoordTransform ) THEN
        CALL CoordinateTransformation( Mesh, CoordTransform, Params, &
            IrreversibleTransformation = .FALSE. )
      END IF

      CALL Info('PartitionMesh','Using partitioning method: '//TRIM(SetMethod))
      
      PartOffset = MAXVAL( ElementPart )       
      CALL Info('PartitionMesh','Partitioning offset: '//TRIM(I2S(PartOffset)))

      SELECT CASE( SetMethod ) 
        
        !CASE( 'metis recursive' ) 
        !CASE( 'metis kway' ) 
        !CASE( 'metis nodal' ) 
        !CASE( 'metis dual' ) 

        CASE( 'directional')
          IF( SetNodes ) THEN
            CALL ClusterNodesByDirection(Params,&
                Mesh,ElementPart,PartitionCand)
          ELSE
            CALL ClusterElementsByDirection(Params,&
                Mesh,ElementPart,PartitionCand)
          END IF

        CASE( 'uniform' )          
          CALL ClusterElementsUniform(Params,&
              Mesh,ElementPart,PartitionCand)
          
        CASE DEFAULT
          CALL Fatal('PartitionMesh','Unspecificed partitioning: '//TRIM(SetMethod))
          
      END SELECT

      IF( PartOffset > 0 ) THEN
        WHERE( PartitionCand ) ElementPart = ElementPart + PartOffset
      END IF


      CALL Info('PartitionMesh','Partitioning of set finished')

      IF( GotCoordTransform ) THEN
        CALL BackCoordinateTransformation( Mesh, DeleteTemporalMesh = .TRUE. )
      END IF

  
    END SUBROUTINE PartitionMeshPart
      

    ! Given a partitioning create a list of Neighbours needed for the communication
    !------------------------------------------------------------------------------

    SUBROUTINE CreateNeighbourList()

      INTEGER :: i,j,k,l,n,m,Partition,lsum,lmax
      INTEGER :: TmpNeighbours(100)
      TYPE(Element_t), POINTER :: Element
      INTEGER :: allocstat

      CALL Info('PartitionMesh','Creating neighbour list for parallel saving')

      n = Mesh % NumberOfNodes
      ALLOCATE( NeighbourList(n) , STAT=allocstat ) 
      IF( allocstat /= 0 ) THEN
        CALL Fatal('PartitionMesh','Allocation error for NeighbourList')
      END IF


      DO i=1,n
        NULLIFY( NeighbourList(i) % Neighbours )
      END DO
      
      DO i=1,Mesh % NumberOfBulkElements
        Element => Mesh % Elements(i)
        Partition = ElementPart(i)
        m = Element % TYPE % NumberOfNodes
        DO j=1,m
          k = Element % NodeIndexes(j)
          IF( .NOT. ASSOCIATED( NeighbourList(k) % Neighbours ) ) THEN
            ALLOCATE( NeighbourList(k) % Neighbours(1), STAT = allocstat )
            IF( allocstat /= 0 ) THEN
              CALL Fatal('PartitionMesh','Allocation error for Neighbours')
            END IF            
            NeighbourList(k) % Neighbours(1) = Partition
          ELSE IF( .NOT. ANY( NeighbourList(k) % Neighbours == Partition ) ) THEN
            l = SIZE( NeighbourList(k) % Neighbours )

            TmpNeighbours(1:l) = NeighbourList(k) % Neighbours(1:l)
            DEALLOCATE( NeighbourList(k) % Neighbours )

            ALLOCATE( NeighbourList(k) % Neighbours(l+1), STAT = allocstat )
            IF( allocstat /= 0 ) THEN
              CALL Fatal('PartitionMesh','Allocation error for Neighbours')
            END IF                       
            NeighbourList(k) % Neighbours(1:l) = TmpNeighbours(1:l)
            NeighbourList(k) % Neighbours(l+1) = Partition
          END IF
        END DO
      END DO

      lmax = 0
      lsum = 0
      DO k=1,n
        l = SIZE(NeighbourList(k) % Neighbours)
        lmax = MAX( lmax, l )
        lsum = lsum + l
      END DO
      
      CALL Info('PartitionMesh','Maximum number of partitions for a node: '//TRIM(I2S(lmax)))
      
      WRITE(Message,'(A,F8.3)') 'Average number of partitiones for a node: ',1.0_dp*lsum/n
      CALL Info('PartitionMesh',Message) 
      
    END SUBROUTINE CreateNeighbourList


  END SUBROUTINE PartitionMeshSolver

