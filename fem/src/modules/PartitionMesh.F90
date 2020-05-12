!------------------------------------------------------------------------------
!> Partition the finite element mesh using a set of different strategies.
!> This is first a solver, but will perhaps eventually be made an internal routine.
!
!> The idea is that one could use recursive strategies for different pieces of the 
!> finite element mesh. This way the physics can be better taken into account than
!> when using stand alone partitioning tools.
!
!> This will hopefully gradually move to the library and become there an fully
!> integrated strategy making this solver obsolete.
!
! P.R.
!------------------------------------------------------------------------------
   SUBROUTINE PartitionMeshSolver( Model,Solver,Timestep,TransientSimulation )
!------------------------------------------------------------------------------
     USE DefUtils
     USE MeshPartition

     IMPLICIT NONE
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver          !< Linear & nonlinear equation solver options
     TYPE(Model_t), TARGET :: Model    !< All model information (mesh, materials, BCs, etc...)
     REAL(KIND=dp) :: Timestep         !< Timestep size for time dependent simulations
     LOGICAL :: TransientSimulation    !< Steady state or transient simulation
!------------------------------------------------------------------------------
     TYPE(ValueList_t), POINTER :: Params
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER, POINTER :: ElementPart(:)
     INTEGER :: n
     LOGICAL :: Found
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


     CALL PartitionMeshSerial( Model, Mesh, Params )
     ElementPart => Mesh % RePartition
     
     CALL CreateNeighbourList()
     
     DirectoryName = ListGetString( Params,'Output Directory',Found )
     IF( Found ) THEN
       CALL WriteMeshToDiskPartitioned( Model, Mesh,DirectoryName, &
           ElementPart, NeighbourList )
     END IF

     IF( ListGetLogical( Params,'Partitioning Exit',Found ) ) STOP

     CALL Info('PartitionMesh','Create partitioning variable')
     CALL SetPartitionVariable()

     CALL Info('PartitionMesh','All done for now')


   CONTAINS


     ! Create a variable for the output of the partitioning.
     ! This is an elemental field, not a nodal one. 
     !------------------------------------------------------
     SUBROUTINE SetPartitionVariable()
       
       TYPE(Variable_t), POINTER :: Var
       INTEGER :: t
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
      
      WRITE(Message,'(A,F8.3)') 'Average number of partitions for a node: ',1.0_dp*lsum/n
      CALL Info('PartitionMesh',Message) 
      
    END SUBROUTINE CreateNeighbourList


  END SUBROUTINE PartitionMeshSolver

