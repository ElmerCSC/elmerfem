!/*****************************************************************************/
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

!> \ingroup ElmerLib 
!> \{

!------------------------------------------------------------------------------
!> Module containing various clustering methods that may be used in conjunction
!> with algebraic multigrid methods or mesh partitioning, for example.
! author Peter Råback
!------------------------------------------------------------------------------


MODULE ClusteringMethods
  
  USE Lists
!  USE CRSMatrix
  USE ElementUtils, ONLY : TangentDirections
  
  IMPLICIT NONE

CONTAINS


!---------------------------------------------------------------
!> This partitions the mesh into a given number of partitions in each 
!> direction. It may be used in clustering multigrid or similar, 
!> and also to internal partitioning within ElmerSolver. 
!---------------------------------------------------------------
  SUBROUTINE ClusterNodesByDirection(Params,Mesh,Clustering,MaskActive)
 
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
    CHARACTER(*), PARAMETER :: Caller="ClusterNodesByDirection"
!---------------------------------------------------------------

    CALL Info(Caller,'Starting clustering',Level=30)

    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = COUNT( MaskActive )
    ELSE
      nsize = Mesh % NumberOfNodes
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal(Caller,'No parameter list associated')
    END IF

    dim = Mesh % MeshDim
    Parray => ListGetConstRealArray( Params,'Clustering Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal(1) = 1.0
      Normal(2) = 1.0d-2
      IF( dim == 3) Normal(3) = 1.0d-4
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
          Divisions(1) = ( nsize / clustersize ) ** 0.5_dp
          Divisions(2) = ( nsize / ( clustersize * Divisions(1) ) )
        ELSE
          Divisions(1:2) = ( nsize / clustersize ) ** (1.0_dp / 3 )
          Divisions(3) = ( nsize / ( clustersize * Divisions(1) * Divisions(2) ) )
        END IF
      ELSE
        CALL Fatal(Caller,'Clustering Divisions not given!')
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
    CALL Info(Caller,Message)
    WRITE(Message,'(A,T25,I10)') 'Max nodes in cluster:',maxpart
    CALL Info(Caller,Message)
    WRITE(Message,'(A,T28,F10.2)') 'Average nodes in cluster:',avepart
    CALL Info(Caller,Message)
    WRITE(Message,'(A,T28,F10.2)') 'Deviation of nodes:',devpart
    CALL Info(Caller,Message)
    

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
        avepart,devpart, dist
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    CHARACTER(*), PARAMETER :: Caller="ClusterNodesByDirection"
!---------------------------------------------------------------

    CALL Info(Caller,'Starting Clustering',Level=30)

    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = COUNT( MaskActive ) 
    ELSE
      nsize = Mesh % NumberOfBulkElements
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal(Caller,'No parameter list associated')
    END IF

    dim = Mesh % MeshDim
    Parray => ListGetConstRealArray( Params,'Clustering Normal Vector',GotIt )
    IF(.NOT. GotIt) THEN
      Parray => ListGetConstRealArray( Params,'Partitioning Normal Vector',GotIt )
    END IF
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal(1) = 1.0
      Normal(2) = 1.0d-2
      IF( dim == 3) THEN
        Normal(3) = 1.0d-4
      ELSE
        Normal(3) = 0.0_dp
      END IF
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
          Divisions(1) = ( nsize / clustersize ) ** 0.5_dp
          Divisions(2) = ( nsize / ( clustersize * Divisions(1) ) )
        ELSE
          Divisions(1:2) = ( nsize / clustersize ) ** (1.0_dp / 3 )
          Divisions(3) = ( nsize / ( clustersize * Divisions(1) * Divisions(2) ) )
        END IF
      ELSE
        CALL Fatal(Caller,'Clustering Divisions not given!')
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

      ! Now compute the weights for each node
      !----------------------------------------
      j = 0
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

        j = j + 1
        Arrange(j) = SUM( Weights * Coord )

        ! Initialize ordering for the current direction
        Order(j) = j
      END DO

      ! Order the distances for given direction, only the active ones
      !--------------------------------------------------------------
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

        ! Now set the partition 
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
    CALL Info(Caller,Message)
    WRITE(Message,'(A,T25,I10)') 'Max nodes in cluster:',maxpart
    CALL Info(Caller,Message)
    WRITE(Message,'(A,T28,F10.2)') 'Average nodes in cluster:',avepart
    CALL Info(Caller,Message)
    WRITE(Message,'(A,T28,F10.2)') 'Deviation of nodes:',devpart
    CALL Info(Caller,Message)
    
    
    IF( ASSOCIATED(Clustering)) THEN
      IF( PRESENT( MaskActive ) ) THEN
        j = 0
        DO i=1, SIZE(MaskActive)
          IF( MaskActive(i) ) THEN
            j = j + 1
            Clustering(i) = Nodepart(j)
          END IF
        END DO
      ELSE
        Clustering = Nodepart 
      END IF
      DEALLOCATE(Nodepart)
    ELSE
      Clustering => Nodepart
      NULLIFY( Nodepart ) 
    END IF
    
    DEALLOCATE(Order,Arrange,NoPart)

  END SUBROUTINE ClusterElementsByDirection



  SUBROUTINE ClusterElementsUniform(Params,Mesh,Clustering,MaskActive,PartitionDivisions)
 
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Clustering(:)
    LOGICAL, OPTIONAL :: MaskActive(:)
    INTEGER, OPTIONAL :: PartitionDivisions(3)
!---------------------------------------------------------------
    LOGICAL :: MaskExists,UseMaskedBoundingBox,Found
    INTEGER :: i,j,k,ind,n,dim,nsize,nmask,clusters
    INTEGER, POINTER :: Iarray(:),ElemPart(:)
    INTEGER, ALLOCATABLE :: NoPart(:)
    INTEGER :: Divisions(3),minpart,maxpart,Inds(3)
    REAL(KIND=dp) :: Coord(3), Weights(3), avepart,devpart
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: BoundingBox(6)
    INTEGER, ALLOCATABLE :: CellCount(:,:,:)
    LOGICAL, ALLOCATABLE :: NodeMask(:)
    CHARACTER(*),PARAMETER :: Caller="ClusterElementsUniform"

    CALL Info(Caller,'Clustering elements uniformly in bounding box',Level=12)

    IF( Mesh % NumberOfBulkElements == 0 ) RETURN
    
    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = SIZE( MaskActive ) 
      nmask = COUNT( MaskActive ) 
      CALL Info(Caller,'Applying division to masked element: '//I2S(nmask),Level=8)
    ELSE
      nsize = Mesh % NumberOfBulkElements 
      nmask = nsize
      CALL Info(Caller,'Applying division to all bulk elements: '//I2S(nsize),Level=8)
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal(Caller,'No parameter list associated')
    END IF

    dim = Mesh % MeshDim

    ! We can use the masked bounding box
    UseMaskedBoundingBox = .FALSE.
    IF( MaskExists ) UseMaskedBoundingBox = ListGetLogical( Params,&
        'Partition Masked Bounding Box',Found ) 

    IF( UseMaskedBoundingBox ) THEN
      ALLOCATE( NodeMask( Mesh % NumberOfNodes ) )
      NodeMask = .FALSE.

      ! Add all active nodes to the mask
      DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        IF( .NOT. MaskActive( i ) ) CYCLE        
        Element => Mesh % Elements(i)
        NodeIndexes => Element % NodeIndexes 
        NodeMask( NodeIndexes ) = .TRUE.
      END DO

      i = COUNT( NodeMask ) 
      CALL Info(Caller,'Masked elements include nodes: '//I2S(i),Level=8)
      
      ! Define the masked bounding box
      BoundingBox(1) = MINVAL( Mesh % Nodes % x, NodeMask )
      BoundingBox(2) = MAXVAL( Mesh % Nodes % x, NodeMask )
      BoundingBox(3) = MINVAL( Mesh % Nodes % y, NodeMask )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % y, NodeMask )
      BoundingBox(5) = MINVAL( Mesh % Nodes % z, NodeMask )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z, NodeMask )

      DEALLOCATE( NodeMask ) 
    ELSE      
      BoundingBox(1) = MINVAL( Mesh % Nodes % x )
      BoundingBox(2) = MAXVAL( Mesh % Nodes % x )
      BoundingBox(3) = MINVAL( Mesh % Nodes % y )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % y )
      BoundingBox(5) = MINVAL( Mesh % Nodes % z )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
    END IF
      
    
    IF( PRESENT( PartitionDivisions ) ) THEN
      Divisions = PartitionDivisions
    ELSE      
      Iarray => ListGetIntegerArray( Params,'Partitioning Divisions',Found)
      IF(.NOT. Found ) THEN
        CALL Fatal(Caller,'> Partitioning Divisions < not given!')
      END IF      
      Divisions = 1
      IF( Found ) THEN
        n = MIN( SIZE(Iarray), dim ) 
        Divisions(1:n) = Iarray(1:n)
      END IF
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

    ALLOCATE(ElemPart(nsize),NoPart(Clusters))
    NoPart = 0
    ElemPart = 0

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
      
      ! Find the center of the element
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
      ElemPart(i) = ind
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

    CALL Info(Caller,'Number of partitions: '//I2S(n),Level=8)
    CALL Info(Caller,'Min elements in cluster: '//I2S(minpart),Level=8)
    CALL Info(Caller,'Max elements in cluster: '//I2S(maxpart),Level=8)

    WRITE(Message,'(A,F10.2)') 'Average elements in cluster:',avepart
    CALL Info(Caller,Message,Level=8)    
    WRITE(Message,'(A,F10.2)') 'Average deviation in size:',devpart
    CALL Info(Caller,Message,Level=8)

    ! Renumber the partitions using only the active ones
    n = 0
    DO i=1,clusters
      IF( NoPart(i) > 0 ) THEN
        n = n + 1
        NoPart(i) = n
      END IF
    END DO
    
    ! Renumbering only needed if there are empty cells
    IF( n < clusters ) THEN
      DO i=1,nsize
        j = ElemPart(i)
        IF( j > 0 ) ElemPart(i) = NoPart(j) 
      END DO
    END IF

    !DO i=1,clusters
    !  PRINT *,'count in part:',i,COUNT( ElemPart(1:nsize) == i ) 
    !END DO
    
    IF( ASSOCIATED( Clustering ) ) THEN
      WHERE( ElemPart > 0 ) Clustering = ElemPart
      DEALLOCATE( ElemPart ) 
    ELSE
      Clustering => ElemPart
      NULLIFY( ElemPart ) 
    END IF
    
    DEALLOCATE(NoPart,CellCount)

    CALL Info(Caller,'Clustering of elements finished',Level=10)

  END SUBROUTINE ClusterElementsUniform

  
!------------------------------------------------------------------------------
!> Subroutine creates clusters of connections in different ways. The default method
!> uses only information from the matrix connectios. The geometric version makes
!> clusters using recursive splitting in each coordinate direction. There is also
!> a version that assumes that the initial mesh was created using extrusion and 
!> this extrusion may be taken back in clustering.
!------------------------------------------------------------------------------
  SUBROUTINE ChooseClusterNodes(Amat, Solver, Components, EliminateDir, CF)

    USE MeshUtils    
    
    TYPE(Matrix_t), POINTER  :: Amat
    TYPE(solver_t), TARGET :: Solver
    INTEGER :: Components
    LOGICAL :: EliminateDir
    INTEGER, POINTER :: CF(:)

    INTEGER :: matsize, nodesize, Component1
    LOGICAL, POINTER :: Bonds(:), Fixed(:),Passive(:)
    INTEGER, POINTER :: Cols(:),Rows(:),CFLayer(:)
    LOGICAL :: GotIt

    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: MaskPerm(:)
    LOGICAL, POINTER :: MaskActive(:)
    CHARACTER(:), ALLOCATABLE :: ClusterMethod 

    ClusterMethod = ListGetString( Solver % Values,'MG Cluster Method',GotIt)
    IF(.NOT. GotIt ) ClusterMethod = 'default'

    SELECT CASE ( ClusterMethod ) 
      
    CASE( 'geometric' )      
      CALL Info('ChooseClusterNodes','Using geometric clustering')
      PSolver => Solver
      Mesh => Solver % Mesh 
      MaskPerm => Solver % Variable % Perm
      IF( ASSOCIATED( MaskPerm ) ) THEN
        ALLOCATE( MaskActive( SIZE( MaskPerm ) ) ) 
        MaskActive = ( MaskPerm > 0 )
      END IF

      IF( Amat % NumberOfRows /= Components * Mesh % NumberOfNodes ) THEN
        CALL Fatal('ChoosClusterNodes','Mismatch in dimensions, geometric clustering works only for two levels!')
      END IF

      CALL ClusterNodesByDirection(Solver % Values,Mesh,CF,MaskActive)

      IF( ASSOCIATED( MaskPerm ) ) DEALLOCATE( MaskActive ) 
    
    CASE( 'unextrude' )
      CALL Info('ChooseClusterNodes','Using dimensional reduction of extruded meshes for clustering')
      Mesh => Solver % Mesh 
      MaskPerm => Solver % Variable % Perm
      IF( Amat % NumberOfRows /= Components * Mesh % NumberOfNodes ) THEN
        PRINT *,'sizes:',Amat % NumberOfRows, Components, Mesh % NumberOfNodes
        CALL Fatal('ChoosClusterNodes','Mismatch in dimensions, extruded node clustering works only for two levels!')
      END IF
      CALL ClusterExtrudedMesh()

    CASE( 'edge' )
      CALL Info('ChooseClusterNodes','Using dimensional edge reduction of extruded edge meshes')
      Mesh => Solver % Mesh 
      MaskPerm => Solver % Variable % Perm
      
      IF( Amat % NumberOfRows /= Components * Mesh % NumberOfEdges ) THEN
        PRINT *,'sizes:',Amat % NumberOfRows, Components, Mesh % NumberOfEdges
        CALL Fatal('ChoosClusterNodes','Mismatch in dimensions, extruded edge clustering works only for two levels!')
      END IF

      CALL ClusterExtrudedEdges()

      
    CASE( 'hybrid' )
      CALL Info('ChooseClusterNodes','Using hybrid clustering method')       
      Component1 = 1
      IF(Components > 1) THEN
        Component1 = ListGetInteger(Solver % Values,'MG Determining Component',&
            GotIt,minv=1,maxv=Components)
        IF(.NOT. GotIt) Component1 = 1
      END IF
      
      Rows => Amat % Rows
      Cols => Amat % Cols
      matsize = Amat % NumberOfRows
      nodesize = matsize / Components
      ALLOCATE( Bonds(SIZE(Amat % Cols)), Passive(nodesize),Fixed(nodesize), CF(nodesize) )

      PSolver => Solver
      Mesh => Solver % Mesh 
      MaskPerm => Solver % Variable % Perm
      IF( Amat % NumberOfRows /= Components * Mesh % NumberOfNodes ) THEN
        CALL Fatal('ChoosClusterNodes','Mismatch in dimensions, geometric clustering works only for two levels!')
      END IF
       
      CALL PassiveExtrudedMesh()

      CALL CMGBonds(Amat, Bonds, Passive, Fixed, Components,Component1)      
      CALL CMGClusterForm(Amat, Bonds, Passive, Fixed, Components, Component1, CFLayer)
      DEALLOCATE(Bonds, Fixed)

      CALL Info('ChooseClusterNodes','Using dimensional reduction of extruded meshes for clustering')
      Mesh => Solver % Mesh 
      MaskPerm => Solver % Variable % Perm
      IF( Amat % NumberOfRows /= Components * Mesh % NumberOfNodes ) THEN
        CALL Fatal('ChoosClusterNodes','Mismatch in dimensions, extruded clustering works only for two levels!')
      END IF
      CALL ClusterExtrudedMesh( CFLayer)

    CASE DEFAULT
      CALL Info('ChooseClusterNodes','Using clustering based on matrix connections')       
      Component1 = 1
      IF(Components > 1) THEN
        Component1 = ListGetInteger(Solver % Values,'MG Determining Component',&
            GotIt,minv=1,maxv=Components)
        IF(.NOT. GotIt) Component1 = 1
      END IF
      
      Rows => Amat % Rows
      Cols => Amat % Cols
      matsize = Amat % NumberOfRows
      nodesize = matsize / Components
      ALLOCATE( Bonds(SIZE(Amat % Cols)), Passive(nodesize),Fixed(nodesize), CF(nodesize) )

      ! For parallel cases we cluster only the "own" dofs and communicate the clustering
      ! afterwards. The not own dofs are skipped.
      Passive = .FALSE.
      CALL SetParallelPassive()

      CALL CMGBonds(Amat, Bonds, Passive, Fixed, Components,Component1)      

      CALL CMGClusterForm(Amat, Bonds, Passive, Fixed, Components, Component1, CF)

      DEALLOCATE(Bonds, Fixed)

    END SELECT

    IF(.FALSE.) CALL Info('ChooseClusterNodes','Clusters chosen')

  CONTAINS


    !-------------------------------------------------------------------------------
    !> This subroutine sets nodes that are now owned by the partition to be passive
    !> before making the clustering. The idea is to communicate the clustering on 
    !> parallel level based on the ownership.
    !-------------------------------------------------------------------------------
    SUBROUTINE SetParallelPassive()
      INTEGER :: i,j

      ! If not parallel then return
      IF( ParEnv % PEs <= 1 ) RETURN

      ! Note that the ParallelInfo for matrices equals the size of the matrix
      ! while for Meshes in equals the number of nodes. Here the former is used.
      DO i=1,nodesize
        j = Components*(i-1) + Component1
        IF( Amat % ParallelInfo % GInterface(j) ) CYCLE        
        IF( Amat % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % Mype ) &
            Passive(i) = .TRUE.
      END DO     
    END SUBROUTINE SetParallelPassive



    !-------------------------------------------------------------------------------
    !> Sets a marker for passive dofs everywhere else than at the top layer.
    !> This is intended for the hybrid method where a cluster is first made on
    !> a level and this is then inherited on the extruded levels. 
    !------------------------------------------------------------------------
    SUBROUTINE PassiveExtrudedMesh()

      INTEGER, POINTER :: TopPointer(:)
      INTEGER :: i, j, nsize
      TYPE(Variable_t), POINTER :: ExtVar
      TYPE(Solver_t), POINTER :: Psolver

      PSolver => Solver      
      nsize = SIZE( MaskPerm )

      CALL Info('PassiveExrudedMesh','Setting active nodes on a layer only')

      CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar, &
          TopNodePointer = TopPointer )
      
      DO i=1,nsize
        ! Node at top is active
        IF( TopPointer(i) == i ) CYCLE

        ! All other nodes, if existing, are passive
        ! Here i refers to physical node, and j to corresponding ordering in matrix
        j = MaskPerm( i )
        IF( j == 0 ) CYCLE
        Passive( j ) = .TRUE.
      END DO

      i = COUNT( Passive ) 
      PRINT *,'Number of passive nodes:',i
      PRINT *,'Number of active nodes:',nsize-i

      DEALLOCATE( TopPointer )

    END SUBROUTINE PassiveExtrudedMesh


    !-------------------------------------------------------------------------------
    !> First detects an extruded structure using mesh information, and 
    !> then clusters the nodes along the extruded directions only.
    !-------------------------------------------------------------------------------
    SUBROUTINE ClusterExtrudedMesh( CFLayer )

      INTEGER, POINTER, OPTIONAL :: CFLayer(:)

      INTEGER, POINTER :: TopPointer(:), DownPointer(:)
      INTEGER :: i, j, k, ind, nsize, NoLayers, cluster, NoClusters, ClusterSize
      TYPE(Variable_t), POINTER :: ExtVar
      TYPE(Solver_t), POINTER :: Psolver

      PSolver => Solver
      
      nsize = SIZE( MaskPerm )
      ALLOCATE( CF( nsize) )

      CF = 0

      ClusterSize = ListGetInteger(Solver % Values,'MG Cluster Size',GotIt)
      
      ! Find the extruded structure 
      CALL DetectExtrudedStructure( Mesh, PSolver, ExtVar, &
          TopNodePointer = TopPointer, DownNodePointer = DownPointer, &
          NumberOfLayers = NoLayers )
  
      ! Numbering of the top layer clusters
      IF( PRESENT( CFLayer ) ) THEN
        NoClusters = MAXVAL( CFLayer ) 
        CF = CFLayer
      ELSE
        NoClusters = 0
        DO i=1,nsize
          IF( TopPointer(i) == i ) THEN
            NoClusters = NoClusters + 1
            CF( MaskPerm(i) ) = NoClusters
          END IF
        END DO
      END IF
          
      ! Number the nodes not on the top layer also
      DO i=1,nsize
        IF( TopPointer(i) == 0 ) CYCLE

        ! start from each top layer node
        IF( TopPointer(i) == i ) THEN
          k = 1 
          ind = i

          ! inherit the first cluster from the top layer
          cluster = CF( MaskPerm(ind) )
          DO j = 1,NoLayers
            ind = DownPointer(ind)

            ! when cluster size is full start a new cluster layer with offset of NoClusters
            IF( k == ClusterSize ) THEN
              k = 0
              cluster = cluster + NoClusters
            END IF
            k = k + 1
            CF( MaskPerm(ind) ) = cluster
          END DO
        END IF
      END DO

      DEALLOCATE( DownPointer, TopPointer )

    END SUBROUTINE ClusterExtrudedMesh

    !-------------------------------------------------------------------------------
    !> First detects an extruded elements using mesh information, and 
    !> then clusters the edges to those of bottom elements,
    !-------------------------------------------------------------------------------
    SUBROUTINE ClusterExtrudedEdges( CFLayer )

      INTEGER, POINTER, OPTIONAL :: CFLayer(:)

      INTEGER, POINTER :: TopPointer(:), BotPointer(:)
      INTEGER :: i, j, k, j2, nsize, NoLayers, NoClusters
      TYPE(Variable_t), POINTER :: ExtVar
      TYPE(Solver_t), POINTER :: Psolver
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Element_t), POINTER :: Element, Element2
      INTEGER, POINTER :: Indexes(:), Indexes2(:), Perm(:), CFPerm(:)
      INTEGER :: Indx, elem, t, t2, dofs, n, n2, m, m2, n0
      INTEGER :: Ncount(5)
      
      
      PSolver => Solver
      Mesh => Solver % Mesh 
      Perm => Solver % Variable % Perm

      
      dofs = Solver % Variable % DOFs
      nsize = Solver % Matrix % NumberOfRows / dofs
      n0 = Mesh % NumberOfNodes

      ALLOCATE( CF( nsize), CFPerm(nsize) )
      CF = 0
      CFPerm = 0

      Ncount = 0
      
      ! Find the extruded structure 
      CALL DetectExtrudedElements( Mesh, PSolver, ExtVar, &
          TopElemPointer = TopPointer, BotElemPointer = BotPointer, &
          NumberOfLayers = NoLayers )

      DO elem=1,Solver % NumberOfActiveElements
        
        ! The index of the element to be mapped
        t = Solver % ActiveElements(elem)
        ! The index of the target element
        t2 = BotPointer( t )
        
        ! Don't map indexes to self
        IF( t2 == t ) THEN
          Ncount(1) = Ncount(1) + 1
          CYCLE
        ELSE
          NCount(2) = Ncount(2) +1
        END IF
        
        ! Note that here the treatment of n and Indexes applied only to classical edge elements
        Element => Mesh % Elements( t )
        n = Element % TYPE % NumberOFEdges
        Indexes => Element % EdgeIndexes        
        
        Element2 => Mesh % Elements( t2 )
        n2 = Element2 % TYPE % NumberOFEdges
        Indexes2 => Element2 % EdgeIndexes
        
        m = Element % TYPE % NumberOFNodes
        m2 = Element2 % TYPE % NumberOFNodes

        IF( n /= n2 ) CALL Fatal('ClusterExtrudeEdges','Cannot map different number of dofs!')
        
        ! Here we assume that the numbering of each element is similar!
        DO i=1,n
          j = Perm(Indexes(i)) 
          j2 = Perm(Indexes2(i))

          !PRINT *,'CF:',j,CF(j),j2
          
          IF( CF(j) /= j2 ) THEN
            Ncount(3) = Ncount(3) + 1
            CF(j) = j2
          ELSE
            Ncount(4) = Ncount(4) + 1
          END IF
        END DO

      END DO
        
      CFPerm = 0
      DO i=1,nsize        
        CFPerm( CF(i) ) = 1
      END DO

      k = 0
      DO i=1,nsize
        IF( CFPerm(i) > 0 ) THEN
          k = k + 1
          CFPerm(i) = k
        END IF
      END DO
      CALL Info('ClusterExtrudedEdges',&
          'Number of reduced dofs: '//I2S(k)//' vs. '//I2S(nsize),Level=9)

      WRITE(Message,'(A,F10.3)') 'Coarse dofs reduction factor',1.0_dp *  nsize / k 
      CALL Info('ClusterExtrudedEdges', Message, Level=7)
            
      DO i=1,nsize
        CF(i) = CFPerm(i)
      END DO

      PRINT *,'Ncount:',Ncount(1:4)
      

      DEALLOCATE( BotPointer, TopPointer )

    END SUBROUTINE ClusterExtrudedEdges


!------------------------------------------------------------------------------
!> Mark the strong connections based on the relative absolute magnitude of the 
!> matrix element.
!------------------------------------------------------------------------------
    SUBROUTINE CMGBonds(Amat, Bonds, Passive, Fixed, Components, Component1)
      
      LOGICAL, POINTER :: Bonds(:), Passive(:), Fixed(:)
      TYPE(Matrix_t), POINTER  :: Amat
      INTEGER :: Components, Component1
      
      REAL(KIND=dp) :: StrongLim
      INTEGER :: nods, cnods, matsize, nodesize, maxconn, newbonds, MaxConns, MinConns
      INTEGER :: i,j,k,cj,ci,no,elimnods,strongbonds,nodeind,matind
      INTEGER, POINTER :: Cols(:),Rows(:),localind(:)
      LOGICAL :: ElimDir, Choose, GotIt
      REAL(KIND=dp), POINTER :: Values(:), measures(:)
      REAL(KIND=dp) :: maxbond, diagbond, dirlim, meas
      
      StrongLim = ListGetConstReal(Solver % Values,'MG Strong Connection Limit',GotIt)
      IF(.NOT. GotIt) StrongLim = 0.06
      MaxConns = ListGetInteger(Solver % Values,'MG Strong Connection Maximum',GotIt)
      MinConns = ListGetInteger(Solver % Values,'MG Strong Connection Minimum',GotIt)
      
      ElimDir = EliminateDir
      DirLim = ListGetConstReal(Solver % Values,'MG Eliminate Dirichlet Limit',GotIt)
      IF(.NOT. GotIt) DirLim = 1.0d-8      
      
      matsize = Amat % NumberOfRows
      Rows   => Amat % Rows
      Cols   => Amat % Cols
      Values => Amat % Values    
      nodesize = matsize / Components
      
      maxconn = 0
      DO matind=1,matsize
        maxconn = MAX(maxconn,Rows(matind+1)-Rows(matind))
      END DO
      ALLOCATE(measures(maxconn),localind(maxconn))
      
      Fixed = .FALSE.
      Bonds = .FALSE.
      strongbonds = 0
      elimnods = 0
      
      DO nodeind=1,nodesize
        
        IF( Passive( nodeind ) ) CYCLE

        matind = (nodeind-1)*Components + Component1
        no = 0
        maxbond = 0.0d0
        
        DO j=Rows(matind),Rows(matind+1)-1
          cj = Cols(j)
          IF( MOD(matind, Components) /= MOD(cj, Components) ) CYCLE
          k = ( cj - 1 )/Components + 1
          IF( Passive( k ) ) CYCLE
          
          IF(cj /= matind ) THEN
            no = no + 1
            localind(no) = j
            measures(no) = ABS( Values(j) )
          ELSE
            diagbond = ABS(Values(j))
          END IF
        END DO
        
        IF(no > 0) maxbond = MAXVAL( measures(1:no) )
        
        ! Mark Dirichlet nodes in order to favor boundaries in future
        IF( ElimDir .AND. maxbond < DirLim * diagbond ) THEN
          Fixed(nodeind) = .TRUE.
          elimnods = elimnods + 1
          CYCLE
        END IF
        
        IF(no == 0) CYCLE
        Choose = .FALSE.
        
        ! Check whether the tentative list would give the right number anaway
        IF(MaxConns > 0 .OR. MinConns > 0) THEN
          newbonds = 0
          DO i=1,no
            IF( measures(i) >= Stronglim * maxbond) newbonds = newbonds + 1
          END DO
          IF(MaxConns > 0 .AND. newbonds > MaxConns) Choose = .TRUE.
          IF(MinConns > 0 .AND. newbonds < MinConns) Choose = .TRUE.
        END IF
        
        IF(.NOT. Choose) THEN
          newbonds = 0
          DO i=1,no
            IF( measures(i) >= StrongLim * maxbond) THEN
              newbonds = newbonds + 1
              Bonds(LocalInd(i)) = .TRUE.
            END IF
          END DO
          strongbonds = strongbonds + newbonds
        ELSE
          CALL SortR(no,LocalInd,measures)
          IF(MaxConns > 0 .AND. newbonds > MaxConns) THEN
            strongbonds = strongbonds + MaxConns
            DO i=1,MaxConns
              Bonds(LocalInd(i)) = .TRUE.
            END DO
          END IF
          
          IF(MinConns > 0 .AND. newbonds < MinConns) THEN
            no = MIN(MinConns, no)
            strongbonds = strongbonds + no
            DO i=1,no
              Bonds(LocalInd(i)) = .TRUE.
            END DO
          END IF
        END IF
      END DO
      
      DEALLOCATE(measures, localind)
      
      IF(elimnods > 0) THEN
        WRITE(Message,'(A,I8)') 'Number of eliminated nodes',elimnods
        CALL Info('CMGBonds',Message)
      END IF
      WRITE(Message,'(A,F8.3)') 'Average number of strong bonds',&
          1.0*strongbonds/nodesize
      CALL Info('CMGBonds',Message)
      
    END SUBROUTINE CMGBonds
    

    !------------------------------------------------------------------------------
    !> Creates clusters of nodes using the given list of strong connections. 
    !> Nodes assigned by the Passive vector may not be included in the 
    !> clusters - others are ignored. The subroutine returns the vector CF which tells
    !> to which cluster each node belongs to.
    !------------------------------------------------------------------------------
    SUBROUTINE CMGClusterForm(Amat, Bonds, Passive, Fixed, Components, Component1, CF)
      
      TYPE(Matrix_t), POINTER :: Amat
      INTEGER :: Components, Component1
      LOGICAL :: Bonds(:),Passive(:),Fixed(:)
      INTEGER, POINTER :: CF(:)
      
      INTEGER :: nods, cnods    
      INTEGER :: i,j,k,l,n,k2,cj,ci,nodeci,nodecj, &
          MaxCon, MaxConInd, RefCon, cind, find, rind, &
          GatMax, anods, points, neworder, oldorder, nodesize,matsize,&
          nextind, sumcon, maxsumcon, ClusterSize, nbonds, &
          nsize, nodeind, nodeind2, matind, BoundPoints, ClusterPoints, &
          ConnectionPoints, ClusterMode, HalfSize
      INTEGER :: LocalCon(100), LocalInd(100)
      REAL(KIND=dp) :: LocalVal(100),LocalVal2(100), maxv
      INTEGER, POINTER :: Con(:)
      INTEGER, POINTER :: Cols(:),Rows(:)
      INTEGER, ALLOCATABLE :: GatLims(:), GatNos(:), ConInd(:), RevConInd(:)
      REAL(KIND=dp), POINTER :: Values(:)
      LOGICAL :: ClusterOrphans, OrphansBest, ClusterGrow, hit, GotIt
      
      
      BoundPoints = ListGetInteger(Solver % Values,'MG Boundary Priority',GotIt)
      IF(.NOT. GotIt) BoundPoints = 1
      ClusterPoints = ListGetInteger(Solver % Values,'MG Cluster Priority',GotIt)
      IF(.NOT. GotIt) ClusterPoints = MIN(1,BoundPoints+1)
      ConnectionPoints = ListGetInteger(Solver % Values,'MG Connection Priority',GotIt)
      IF(.NOT. GotIt) ConnectionPoints = 1
      ClusterSize = ListGetInteger(Solver % Values,'MG Cluster Size',GotIt)
      HalfSize = ListGetInteger(Solver % Values,'MG Cluster Half Size',GotIt)
      IF(.NOT. GotIt) HalfSize = (ClusterSize-1)/2 
      
      ! Particularly for small clusters the orphan control seems to be important
      ClusterOrphans = ListGetLogical(Solver % Values,'MG Cluster Orphans',GotIt)
      IF(.NOT. GotIt) ClusterOrphans = .TRUE.
      OrphansBest = ListGetLogical(Solver % Values,'MG Cluster Orphans Best',GotIt)
      ClusterGrow = ListGetLogical(Solver % Values,'MG Cluster Grow',GotIt)
      
      IF(ClusterSize == 0) THEN
        ClusterMode = 0
      ELSE IF(ClusterSize <= 2) THEN
        ClusterMode = 1
      ELSE 
        ClusterMode = 2
      END IF
      
      matsize = Amat % NumberOfRows
      nodesize = matsize / Components
      Rows   => Amat % Rows
      Cols   => Amat % Cols
      Values => Amat % Values
      
      CF = 0
      ALLOCATE(Con(nodesize))
      Con = 0
      cnods = 0
      cind = 0
      
      ! Set an initial measure of importance
      !--------------------------------------
      DO nodeind = 1, nodesize
        IF( Passive( nodeind ) ) CYCLE
        IF(Fixed(nodeind)) THEN
          Con(nodeind) = 0
        ELSE
          Con(nodeind) = 1
        END IF
      END DO
      
      ! Add a point for each clustered neighbour 
      !-------------------------------------------
      DO nodeind = 1, nodesize
        IF( Passive(nodeind) ) CYCLE
        IF(CF(nodeind) > 0) THEN
          matind = Components*(nodeind-1) + Component1
          DO i=Rows(matind),Rows(matind+1)-1
            IF(.NOT. Bonds(i)) CYCLE
            ci = Cols(i)     
            nodeci = (ci-Component1) / Components + 1
            IF( PASSIVE(nodeci) ) CYCLE
            IF(.NOT. Fixed(nodeci) .AND. CF(nodeci) == 0) THEN
              Con(nodeci) = Con(nodeci) + ClusterPoints
            END IF
          END DO
        END IF
      END DO
      
      ! Add some weight if the node is connected to a Dirichlet node
      ! This way the clustering will follow the natural boundaries
      !-------------------------------------------------------------
      IF(BoundPoints /= 0) THEN
        DO nodeind = 1, nodesize
          IF(Passive(nodeind)) CYCLE
          IF(Fixed(nodeind)) CYCLE
          matind = Components*(nodeind-1) + Component1
          DO i=Rows(matind),Rows(matind+1)-1
            IF(.NOT. Bonds(i)) CYCLE
            ci = Cols(i)
            nodeci = (ci-Component1) / Components + 1
            IF(Passive(nodeci)) CYCLE
            IF(Fixed(nodeci)) THEN
              Con(nodeind) = Con(nodeind) + BoundPoints
            END IF
          END DO
        END DO
      END IF
      
      MaxCon = MAXVAL( Con ) 
      GatMax = 2 * MaxCon + 100
      
      ! Bookkeeping is kept on category boundaries of the points
      !---------------------------------------------------------
      ALLOCATE( GatLims(GatMax), GatNos(GatMax) )
      GatLims = 0
      GatNos = 0
      
      ! Number of points in different categories 
      !-----------------------------------------
      DO nodeind = 1, nodesize
        IF(Passive(nodeind)) CYCLE
        IF(Con(nodeind) == 0) CYCLE
        GatNos(Con(nodeind)) = GatNos(Con(nodeind)) + 1 
      END DO
      anods = SUM(GatNos)
      
      DO nodeind = 2, MaxCon
        GatLims(nodeind) = GatLims(nodeind-1) + GatNos(nodeind-1)
      END DO
      
      ALLOCATE( ConInd(nodesize), RevConInd(anods) )
      ConInd = 0
      RevConInd = 0
      
      GatNos = 0
      DO nodeind = 1, nodesize
        IF(Passive(nodeind)) CYCLE
        IF(Con(nodeind) == 0) CYCLE
        GatNos(Con(nodeind)) = GatNos(Con(nodeind)) + 1
        ConInd(nodeind) = GatNos(Con(nodeind)) + GatLims(Con(nodeind))
      END DO
      
      RevConInd = 0
      DO nodeind = 1, nodesize
        IF(Passive(nodeind)) CYCLE
        IF(Con(nodeind) == 0) CYCLE
        RevConInd(ConInd(nodeind)) = nodeind
      END DO
      
      DO WHILE( MaxCon > 0 ) 
        nodeind = RevConInd(anods)
        
        IF(ClusterSize /= 0) THEN    
          ! Assume that the cluster consists of the point with maximum weight and its
          ! closest tightly coupled neighbours.
          !---------------------------------------------------------------------------
          n = 0
          nbonds = 0
          matind = Components*(nodeind-1)+Component1
          DO j=Rows(matind),Rows(matind+1)-1
            
            IF(Bonds(j)) THEN
              cj = Cols(j)   
              nodecj = (cj-Component1) / Components + 1

              IF( Passive(nodecj) ) CYCLE
              IF( Fixed(nodecj) ) CYCLE

              nbonds = nbonds + 1    
              IF(CF(nodecj) == 0) THEN
                n = n + 1
                LocalCon(n) = Con(nodecj)
                LocalInd(n) = nodecj
                LocalVal(n) = ABS( Values(j) )
              END IF
            END IF
          END DO
          
          ! Favor maximum number of connections before strength of connections
          !-------------------------------------------------------------------
          IF( n > ClusterSize - 1 ) THEN
            
            ! Renormalize the effect of strength of connection to [0,1]
            maxv = MAXVAL(LocalVal(1:n))
            LocalVal(1:n) = 0.999 * ConnectionPoints * LocalVal(1:n) / maxv + LocalCon(1:n)
            k = HalfSize
            IF( HalfSize > 0) THEN
              CALL SortR(n,LocalInd,LocalVal)
            END IF
            
            IF( ClusterMode == 2 ) THEN
              k = HalfSize
              ! Ensure that the strongest connections are always at the top of the list 
              IF(k > 0) LocalVal(1:k) = LocalVal(1:k) + 1000.0d0
              
              ! This and the other commented line could be checked for higher dimension
              IF(.FALSE.) LocalVal(k+1:n) = LocalVal(k+1:n) + 10.0d0
              
              ! Add points for the connection to the already chosen ones
              DO l=1,k
                nodeind2 = LocalInd(l)
                matind = Components*(nodeind2-1)+Component1
                DO j=Rows(matind),Rows(matind+1)-1
                  IF(Bonds(j)) THEN
                    cj = Cols(j)   
                    nodecj = (cj-Component1) / Components + 1
                    IF(Passive(nodecj)) CYCLE

                    ! Let the connection dominate over other points
                    DO k2=k+1,n
                      IF( LocalInd(k2) == nodecj) THEN
                        LocalVal(k2) = LocalVal(k2) + 10.0d0
                        EXIT
                      END IF
                    END DO
                    IF(ClusterGrow .AND. k2 > n) THEN
                      n=n+1
                      LocalInd(n) = nodecj
                      LocalVal(n) = 10.0d0
                    END IF
                  END IF
                END DO
                
              END DO
            END IF
            
            CALL SortR(n,LocalInd,LocalVal)
            n = ClusterSize - 1
          END IF
          

          ! If a tightly coupled node cannot be a start of a new cluster join it 
          ! to the cluster of its most strongly coupled neighbour
          ! The other option is that the orphan node builds its own cluster.
          !----------------------------------------------------------------------
          IF( ClusterOrphans .AND. n == 0 .AND. nbonds > 0) THEN
            
            maxv = 0.0d0
            DO j=Rows(matind),Rows(matind+1)-1
              IF(Bonds(j)) THEN
                cj = Cols(j)
                nodecj = (cj-Component1) / Components + 1               
                IF( Passive(nodecj) ) CYCLE
                IF( Fixed(nodecj) ) CYCLE
                k = CF(nodecj) 
                IF( k > 0) THEN
                  IF( ABS(Values(j)) > maxv) THEN
                    maxv = ABS( Values(j) )
                    rind = k
                  END IF
                END IF
              END IF
            END DO
            
            ! If there are many possible candidate parents for the orphan take into account the 
            ! sum of all normalized contributions
            !----------------------------------------------------------------------------------
            IF( OrphansBest .AND. nbonds > 2) THEN
              LocalInd(1:nbonds) = 0
              LocalVal(1:nbonds) = 0.0_dp 
              
              DO j=Rows(matind),Rows(matind+1)-1
                IF(Bonds(j)) THEN
                  cj = Cols(j)
                  nodecj = (cj-Component1) / Components + 1               

                  IF( Passive(nodecj) ) CYCLE
                  IF( Fixed(nodecj) ) CYCLE
                  k = CF(nodecj) 
                  
                  IF( k > 0) THEN
                    hit = .FALSE.
                    DO k2=1,n              
                      IF( LocalInd(k2) == k) THEN
                        hit = .TRUE.
                        EXIT
                      END IF
                    END DO
                    IF(.NOT. hit) THEN
                      n = n + 1
                      k2 = n
                      LocalInd(k2) = k
                    END IF
                    LocalVal(k2) = LocalVal(k2) + ABS(Values(j)) / maxv + 0.0_dp
                  END IF
                END IF
              END DO
              
              maxv = 0.0_dp
              DO k=1,n
                IF( LocalVal(k) > maxv) THEN
                  maxv = LocalVal(k) 
                  k2 = LocalInd(k)
                END IF
              END DO
              n = 0
              rind = k2
            END IF
          ELSE
            cind = cind + 1
            rind = cind
          END IF
          
          DO j=0,n
            IF(j==0) THEN
              nodecj = nodeind
            ELSE
              nodecj = LocalInd(j)
            END IF
            
            CF(nodecj) = rind
            cnods = cnods + 1
            
            ! Recompute the measure of importance for the neighbours
            cj = Components*(nodecj-1) + Component1
            DO i=Rows(cj),Rows(cj+1)-1
              IF(.NOT. Bonds(i)) CYCLE

              ci = Cols(i)          
              nodeci = (ci-Component1)/Components + 1
              
              IF( Passive( nodeci) ) CYCLE
              IF( Fixed(nodeci) ) CYCLE
              IF( CF(nodeci) /= 0) CYCLE
              
              DO l=1,ClusterPoints
                points = Con(nodeci) + 1
                oldorder = ConInd(nodeci)
                IF(GatLims(points) == 0) GatLims(points) = anods
                neworder = GatLims(points)
                Con(nodeci) = points
                GatLims(points) = GatLims(points) - 1
                
                IF(neworder /= oldorder) THEN
                  k = RevConInd(neworder)
                  
                  ConInd(nodeci) = neworder
                  ConInd(k) = oldorder
                  
                  RevConInd(neworder) = nodeci
                  RevConInd(oldorder) = k
                END IF
              END DO
              
            END DO
          END DO
        END IF  ! ClusterSize /= 0
        
        ! Assume that the center of clustering is the node with maximum weight,
        ! or any of its strongly coupled neighbours. Compute the total weight of
        ! the candidate clusters.       
        ! Compute first the cluster points for the point with maximum points
        !---------------------------------------------------------------------
        IF(ClusterSize == 0) THEN
          sumcon = 0
          nbonds = 0
          matind = Components*(nodeind-1) + Component1
          DO j=Rows(matind),Rows(matind+1)-1
            IF(.NOT. Bonds(j)) CYCLE

            cj = Cols(j)        
            nodecj = (cj-Component1)/Components + 1
            nbonds = nbonds + 1
            IF(CF(nodecj) == 0) sumcon = sumcon + Con(nodecj) 
          END DO
          maxsumcon = sumcon
          nextind = nodeind
          
          IF(ClusterOrphans .AND. sumcon == 0 .AND. nbonds > 0) THEN
            ! If a tightly coupled node cannot be a start of a new cluster join it 
            ! to the cluster of its most strongly coupled neighbour
            maxv = 0.0d0
            DO j=Rows(matind),Rows(matind+1)-1
              IF(.NOT. Bonds(j)) CYCLE
              
              cj = Cols(j)   
              nodecj = (cj-Component1)/Components + 1
              IF( Fixed(nodecj) ) CYCLE
              IF( CF(nodecj) > 0) THEN
                IF( ABS(Values(j)) > maxv) THEN
                  maxv = ABS( Values(j) )
                  rind = CF(nodecj)
                END IF
              END IF
            END DO
          ELSE
            ! Next compute the cluster points for all its strongly coupled 
            ! unclustered neighbours
            nbonds = 0
            DO j=Rows(matind),Rows(matind+1)-1
              IF(.NOT. Bonds(j)) CYCLE

              sumcon = 0
              cj = Cols(j)                  
              nodecj = (cj-Component1)/Components + 1
              IF(.NOT. Fixed(nodecj) .AND. CF(nodecj) == 0) THEN
                DO i=Rows(cj),Rows(cj+1)-1
                  IF(Bonds(i)) THEN
                    ci = Cols(i)          
                    nodeci = (ci-Component1)/Components + 1
                    IF(CF(nodeci) == 0) sumcon = sumcon + Con(nodeci)
                  END IF
                END DO
              END IF
              IF(sumcon > maxsumcon) THEN
                maxsumcon = sumcon
                nextind = nodecj
              END IF
            END DO
            cind = cind + 1
            rind = cind
          END IF
          
          ! The new center for clustering is nextind
          nodeind = nextind
          CF(nodeind) = rind
          cnods = cnods + 1
          
          ! Go through all strongly bonded neighbours to 'ind'
          matind = Components*(nodeind-1) + Component1
          DO j=Rows(matind),Rows(matind+1)-1
            IF(.NOT. Bonds(j)) CYCLE
            cj = Cols(j)
            nodecj = (cj-Component1)/Components + 1
            
            IF( Passive(nodecj) ) CYCLE
            
            IF(CF(nodecj) == 0) THEN
              CF(nodecj) = cind
              cnods = cnods + 1
              IF(Fixed(nodecj)) CYCLE
              
              ! Recompute the measure of importance for the neighbours
              DO i=Rows(cj),Rows(cj+1)-1
                IF(.NOT. Bonds(i)) CYCLE
                ci = Cols(i)          
                nodeci = (ci-Component1)/Components + 1

                IF( Passive(nodeci) ) CYCLE
                IF( Fixed(nodeci) ) CYCLE
                IF( CF(nodeci) /= 0) CYCLE
                
                DO l=1,ClusterPoints
                  points = Con(nodeci) + 1
                  oldorder = ConInd(nodeci)
                  IF(GatLims(points) == 0) GatLims(points) = anods
                  neworder = GatLims(points)
                  Con(nodeci) = points
                  GatLims(points) = GatLims(points) - 1
                  
                  IF(neworder /= oldorder) THEN
                    k = RevConInd(neworder)
                    
                    ConInd(nodeci) = neworder
                    ConInd(k) = oldorder
                    
                    RevConInd(neworder) = nodeci
                    RevConInd(oldorder) = k
                  END IF
                END DO
                  
              END DO
            END IF
          END DO
        END IF  ! ClusterSize == 0
        
        ! Shorten the list from top in case the node is already set
        Refcon = 0
        DO WHILE( anods > 0)
        ! DO WHILE( anods > 0 .AND. CF(RevConInd(anods)) /= 0)
          IF (CF(RevConInd(anods)) == 0) EXIT
          MaxCon = Con(RevConInd(anods))          
          Refcon = MAX(MaxCon,RefCon)
          anods = anods - 1
          Refcon = MAX(MaxCon,RefCon)            
        END DO
        
        IF(anods > 0) THEN
          MaxCon = Con(RevConInd(anods))
          RefCon = MAX(MaxCon,RefCon)
          GatLims((MaxCon+1):(Refcon+1)) = 0
        ELSE        
          MaxCon = 0
        END IF
      END DO
      
      DEALLOCATE(Con, GatLims, GatNos, ConInd, RevConInd)
      
      WRITE(Message,'(A,I8)') 'Number of clusters',cind
      CALL Info('CMGClusterForm',Message)
      WRITE(Message,'(A,F8.3)') 'Average size of clusters',1.0*cnods/cind
      CALL Info('CMGClusterForm',Message)
    
    END SUBROUTINE CMGClusterForm
    
  END SUBROUTINE ChooseClusterNodes



END MODULE ClusteringMethods


!> \} ElmerLib
