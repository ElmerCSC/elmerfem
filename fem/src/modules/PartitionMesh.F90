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
     INTEGER :: i,j,k,n
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
     ALLOCATE( PartitionCand(n), ElementSet(n), ElementPart(n) )
     PartitionCand = .FALSE.
     ElementSet = 0
     ElementPart = 0
     
     n = MAX( Model % NumberOfBCs, Model % NumberOfEquations ) 
     ALLOCATE( ParameterInd(n) ) 
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

     
!     IF( NumberOfBoundarySets > 0 ) THEN
!     CALL Info('PartitionMesh','Extending the boundary patitioning into the mesh') 
!     IF( NumberOfBoundarySets > 0 ) THEN
!       CALL ExtentMeshPart()
!     END IF
     
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

     SUBROUTINE SetPartitionVariable()
       
       TYPE(Variable_t), POINTER :: Var
       INTEGER :: t
       INTEGER, POINTER :: DGIndexes(:)
       TYPE(Element_t), POINTER :: Element
       
       Var => Solver % Variable 
       IF( .NOT. ASSOCIATED( Var ) ) RETURN
       
       CALL Info('PartitionMesh','Setting discontinuous field > Partition < ')
       
       Var % Values = 0.0_dp
       
       DO t=1,Mesh % NumberOfBulkElements 
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
       
       
       SeparateBoundarySets = ListGetLogical( Params, &
           'Partitioning Separate Boundary Set', Found)
       
       ALLOCATE( BCPart( Model % NUmberOfBCs ) ) 
       BCPart = 0
RETURN

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
             ElementSet( i ) = BCPart( bc_id )
             EXIT
           END IF
         END DO
       END DO

       CALL Info('PartitionMesh','Number of sets for boundary partitioning: '//TRIM(I2S(NumberOfParts)))
      
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
                    
       ElementSet = 0 
       
       ALLOCATE( EquationPart( Model % NumberOfEquations ) ) 
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
      
      
      PartitionCand = ( ElementSet == SetNo )
      n = Solver % Mesh % NumberOfBulkElements
      
      NoCandElements = COUNT( PartitionCand ) 

      CALL Info('PartitionMesh','Number of elements in set: '//TRIM(I2S(NoCandElements)))

      IF( NoCandElements == 0 ) RETURN
         
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

      CALL Info('PartitionMesh','Partitioning of set finished')

      IF( GotCoordTransform ) THEN
        CALL BackCoordinateTransformation( Mesh, DeleteTemporalMesh = .TRUE. )
      END IF
      
  
    END SUBROUTINE PartitionMeshPart
      

    SUBROUTINE CreateNeighbourList()

      INTEGER :: i,j,k,l,n,m,Partition,lsum,lmax
      INTEGER :: TmpNeighbours(100)
      TYPE(Element_t), POINTER :: Element

      CALL Info('PartitionMesh','Creating neighbour list for parallel saving')

      n = Mesh % NumberOfNodes
      ALLOCATE( NeighbourList(n) ) 

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
            ALLOCATE( NeighbourList(k) % Neighbours(1) ) 
            NeighbourList(k) % Neighbours(1) = Partition
          ELSE IF( .NOT. ANY( NeighbourList(k) % Neighbours == Partition ) ) THEN
            l = SIZE( NeighbourList(k) % Neighbours )

            TmpNeighbours(1:l) = NeighbourList(k) % Neighbours(1:l)
            DEALLOCATE( NeighbourList(k) % Neighbours )

            ALLOCATE( NeighbourList(k) % Neighbours(l+1) ) 
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

