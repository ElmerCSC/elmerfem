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
! *  Authors: Peter Råback & Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi & Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 16.6.2011
! *
! *****************************************************************************/



!> \ingroup ElmerLib
!> \{

!-----------------------------------------------------------------------------
!> Module including generic utilities for particle dynamics and tracking.
!-----------------------------------------------------------------------------
 
MODULE ParticleUtils
  
  USE DefUtils
  USE Lists
  USE MeshUtils
  USE GeneralUtils
  USE SaveUtils
  
  IMPLICIT NONE

  TYPE Particle_t
    INTEGER :: Dim, NumberOfParticles=0, MaxNumberOfParticles=0, &
               NumberOfMovingParticles = 0, &
               TimeOrder = 0, FirstGhost = 0, NumberOfGroups = 0
    TYPE(Variable_t), POINTER :: Variables => NULL()	
    
    REAL(KIND=dp) :: time, dtime
    LOGICAL :: DtConstant = .TRUE., RK2 = .FALSE.  
    INTEGER :: DtSign = 1    

    REAL(KIND=dp), POINTER :: Coordinate(:,:) => NULL()
    REAL(KIND=dp), POINTER :: PrevCoordinate(:,:) => NULL()
    REAL(KIND=dp), POINTER :: Velocity(:,:) => NULL()
    REAL(KIND=dp), POINTER :: PrevVelocity(:,:) => NULL()
    REAL(KIND=dp), POINTER :: Force(:,:) => NULL()
    REAL(KIND=dp), POINTER :: uvw(:,:) => NULL()
    
    INTEGER, POINTER :: FaceIndex(:) => NULL()
    INTEGER, POINTER :: Status(:) => NULL()
    INTEGER, POINTER :: ElementIndex(:) => NULL()
    INTEGER, POINTER :: NodeIndex(:) => NULL()
    INTEGER, POINTER :: Partition(:) => NULL()
    INTEGER, POINTER :: Group(:) => NULL()
    
    ! Data structure for the particle-particle interaction 
    INTEGER :: MaxClosestParticles
    LOGICAL :: NeighbourTable = .FALSE.   
    INTEGER, POINTER :: ClosestNode(:) => NULL()
    INTEGER, POINTER :: ClosestParticle(:) => NULL()
    INTEGER, POINTER :: NoClosestParticle(:) => NULL()
    INTEGER, POINTER :: CumClosestParticle(:) => NULL()
    
    ! Mark the internal elements without any interface nodes
    LOGICAL, POINTER :: InternalElements(:) => NULL()
    
    ! Local and global bounding boxes
    REAL(KIND=dp) :: LocalMinCoord(3), LocalMaxCoord(3)
    REAL(KIND=dp) :: GlobalMinCoord(3), GlobalMaxCoord(3)
    
  END TYPE Particle_t
 

  INTEGER, PARAMETER :: &
      PARTICLE_ALLOCATED = 1, &
      PARTICLE_WAITING = 2, &
      PARTICLE_INITIATED = 3, &
      PARTICLE_LOCATED = 4, &
      PARTICLE_MOVING = 5, & 
      PARTICLE_FACEBOUNDARY = 6, &
      PARTICLE_PARTBOUNDARY = 7, &
      PARTICLE_HIT = 8, & 
      PARTICLE_READY = 9, &
      PARTICLE_FIXEDCOORD = 10, &
      PARTICLE_FIXEDVELO = 11, &
      PARTICLE_WALLBOUNDARY = 12, &
      PARTICLE_LOST = 13, &
      PARTICLE_GHOST = 14


  TYPE(Particle_t), TARGET, SAVE :: GlobalParticles
 
  

CONTAINS


!----------------------------------------------------------------------------
!> Counts particles in different categories.
!----------------------------------------------------------------------------
  SUBROUTINE ParticleStatusCount(Particles)
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: i,j,k,NoParticles
    INTEGER :: StatusCount(PARTICLE_GHOST)
    CHARACTER (len=12), PARAMETER :: StatusString(14) = [ &
        "Allocated   ", &
 	"Waiting     ", &
	"Initiated   ", &
	"Located     ", &
	"Moving      ", &
	"FaceBoundary", &
        "PartBoundary", &
	"Hit         ", &
	"Ready       ", &
	"FixedCoord  ", &
	"FixedVelo   ", &
	"WallBoundary", &
        "Lost        ", &
	"Ghost       "]

    StatusCount = 0
    NoParticles = Particles % NumberOfParticles
    DO i=1,NoParticles
      j = Particles % Status(i)
      StatusCount( j ) = StatusCount( j ) + 1
    END DO
      
    CALL Info('ParticleStatusCount','Information on particle status:')
    k = ParallelReduction( NoParticles ) 
    WRITE(Message,'(A,T18,I0)') 'Total: ',k
    CALL Info('ParticleStatusCount',Message,Level=8)
    DO i=1,PARTICLE_GHOST
      j = StatusCount(i)
      k = ParallelReduction( j ) 
      IF( k == 0 ) CYCLE
      WRITE(Message,'(A,T18,I0)') TRIM(StatusString(i))//': ',k
      CALL Info('ParticleStatusCount',Message,Level=8)
    END DO
    
  END SUBROUTINE ParticleStatusCount



  !---------------------------------------------------------
  ! The following subroutines make the data structure 
  ! transparent in the user subrouines and thereby make
  ! them more recilient to time.
  !> Returns coordinates of the particle.
  !---------------------------------------------------------
  FUNCTION GetParticleCoord(Particles,No) RESULT ( Coord )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Coord(3) 
    
    INTEGER :: dim
    
    Coord(3) = 0.0_dp
    dim = Particles % dim
    Coord(1:dim) = Particles % Coordinate(no,1:dim)
  END FUNCTION GetParticleCoord

  !> Returns velocity of the particle.
  !---------------------------------------------------------
  FUNCTION GetParticleVelo(Particles,No) RESULT ( Coord )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Coord(3) 
    
    INTEGER :: dim
    
    Coord(3) = 0.0_dp
    dim = Particles % dim
    Coord(1:dim) = Particles % Velocity(no,1:dim)
  END FUNCTION GetParticleVelo

  !> Returns force acting on the particle.
  !-------------------------------------------------------
  FUNCTION GetParticleForce(Particles,No) RESULT ( Coord )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Coord(3) 
    
    INTEGER :: dim
    
    Coord(3) = 0.0_dp
    dim = Particles % dim
    Coord(1:dim) = Particles % Force(no,1:dim)
  END FUNCTION GetParticleForce

  !> Sets the particle coordinates.
  !--------------------------------------------------------
  SUBROUTINE SetParticleCoord(Particles,No,Coord)
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Coord(3)     
    INTEGER :: dim
    
    dim = Particles % dim    
    Particles % Coordinate(no,1:dim) = Coord(1:dim)
  END SUBROUTINE SetParticleCoord

  !> Sets the particle velocity.
  !-------------------------------------------------------
  SUBROUTINE SetParticleVelo(Particles,No,Velo)
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Velo(3)     
    INTEGER :: dim
    
    dim = Particles % dim    
    Particles % Velocity(no,1:dim) = Velo(1:dim)
  END SUBROUTINE SetParticleVelo

  !> Sets the particle force.
  !-------------------------------------------------------
  SUBROUTINE SetParticleForce(Particles,No,Force)
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Force(3)     
    INTEGER :: dim
    
    dim = Particles % dim    
    Particles % Force(no,1:dim) = Force(1:dim)
  END SUBROUTINE SetParticleForce


  !> Gets the previous particle coordinate.
  !-----------------------------------------------------------
  FUNCTION GetParticlePrevCoord(Particles,No) RESULT ( Coord )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: Coord(3) 
    
    INTEGER :: dim
    
    Coord(3) = 0.0_dp
    dim = Particles % dim

    Coord(1:dim) = Particles % PrevCoordinate(no,1:dim)

  END FUNCTION GetParticlePrevCoord
 
  !> Gets the local coordinates of the element of the given particle.
  !-------------------------------------------------------------------
  SUBROUTINE GetParticleUVW(Particles,No, u, v, w )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: u,v
    REAL(KIND=dp), OPTIONAL :: w
    INTEGER :: dim
    
    dim = Particles % dim
    
    u = Particles % UVW(no,1)
    v = Particles % UVW(no,2)
    IF( PRESENT( w ) ) THEN
      IF( dim == 3 ) THEN
        w = Particles % UVW(no,3)
      ELSE
        w = 0.0_dp
      END IF
    END IF
  END SUBROUTINE GetParticleUVW
  
  !> Sets the local coordinates of the element of the given particle.
  !------------------------------------------------------------------
  SUBROUTINE SetParticleUVW(Particles,No,u,v,w)
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    REAL(KIND=dp) :: u,v
    REAL(KIND=dp), OPTIONAL :: w
    
    INTEGER :: dim
    
    dim = Particles % dim
    
    Particles % UVW(no,1) = u
    Particles % UVW(no,2) = v
    IF( PRESENT( w ) ) THEN
      IF( dim == 3 ) THEN
        Particles % UVW(no,3) = w
      END IF
    END IF
  END SUBROUTINE SetParticleUVW

  !> Adds a displacement to the particle coordinates.
  !------------------------------------------------------------------
  SUBROUTINE AddParticleCoord(Particles,No,Coord)
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No, DerOrder
    REAL(KIND=dp) :: Coord(3)     
    INTEGER :: dim    
    dim = Particles % dim    
    Particles % Coordinate(no,1:dim) = &
        Particles % Coordinate(no,1:dim) + Coord(1:dim)
  END SUBROUTINE AddParticleCoord
  
  !> Adds a velocity difference to the particle velocity.
  !-------------------------------------------------------------------
  SUBROUTINE AddParticleVelo(Particles,No,Coord)
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No, DerOrder
    REAL(KIND=dp) :: Coord(3)     
    INTEGER :: dim    
    dim = Particles % dim
    Particles % Velocity(no,1:dim) = &
        Particles % Velocity(no,1:dim) + Coord(1:dim)
  END SUBROUTINE AddParticleVelo
 
  !> Adds to the force acting on the particle.
  !-------------------------------------------------------------------
   SUBROUTINE AddParticleForce(Particles,No,Force)
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No, DerOrder
    REAL(KIND=dp) :: Force(3)     
    INTEGER :: dim    
    dim = Particles % dim   
    Particles % Force(no,1:dim) = &
        Particles % Force(no,1:dim) + Force(1:dim)
  END SUBROUTINE AddParticleForce

  !> Gets the status of the particle.
  !-------------------------------------------------------------------
  FUNCTION GetParticleStatus(Particles,No) RESULT ( Status )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Status
    
    Status = Particles % Status(No)
  END FUNCTION GetParticleStatus
  
  !> Sets the status of the particle.
  !-------------------------------------------------------------------
  SUBROUTINE SetParticleStatus(Particles,No,Status )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Status
    
    Particles % Status(No) = Status
  END SUBROUTINE SetParticleStatus
  
  !> Gets the elements where the particle is located in, or was located last time.  
  !-------------------------------------------------------------------------------
  FUNCTION GetParticleElement(Particles,No) RESULT ( Index ) 
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Index
    
    Index = Particles % ElementIndex(No)
  END FUNCTION GetParticleElement
  
  !> Sets the element where the particle is located.
  !-------------------------------------------------------------------------------
  SUBROUTINE SetParticleElement(Particles,No,Index )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Index
    
    Particles % ElementIndex(No) = Index
  END SUBROUTINE SetParticleElement
  
  !> Gets the closest node related to the particle.
  !-------------------------------------------------------------------------------
  FUNCTION GetParticleNode(Particles,No) RESULT ( Index ) 
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Index
    
    Index = Particles % NodeIndex(No)
  END FUNCTION GetParticleNode
  
  !> Sets the closest node related to the particle.
  !--------------------------------------------------------------------------------
  SUBROUTINE SetParticleNode(Particles,No,Index )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Index
    
    Particles % NodeIndex(No) = Index
  END SUBROUTINE SetParticleNode
 
  !> Get the group in which the particle belongs to
  !--------------------------------------------------------------------------------
  FUNCTION GetParticleGroup(Particles,No) RESULT ( Index ) 
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Index

    IF( Particles % NumberOfGroups > 0 ) THEN
      Index = Particles % Group(No)
    ELSE
      Index = 0
    END IF
      
  END FUNCTION GetParticleGroup
  
  !> Sets the group in which the particle belongs to 
  !--------------------------------------------------------------------------------
  SUBROUTINE SetParticleGroup(Particles,No,Index )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No
    INTEGER :: Index

    IF( Particles % NumberOfGroups > 0 ) THEN
      Particles % Group(No) = Index
    ELSE
      CALL Warn('SetParticleGroup','Cannot set particle because there is only one group!')
    END IF
    
  END SUBROUTINE SetParticleGroup
 
  !---------------------------------------------------------
  !> The subroutine marks the elements which are not on the 
  !> boundary, either internal or external one. 
  !> This information may be used to speed up different 
  !> loops where particle-boundary interaction is needed.
  !---------------------------------------------------------
  SUBROUTINE MarkInternalElements( Particles )
    
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Element_t), POINTER :: BulkElement, BulkElement2, BoundaryElement
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Body, Body2
    INTEGER :: t,i,j,imax,body_id,body_id2,mat_id,mat_id2,bf_id,bf_id2,dim,istat
    INTEGER :: NumberOfElements
    LOGICAL, POINTER :: InternalElements(:)
    LOGICAL :: Found,Hit
    
    Mesh => GetMesh()
    Dim = Mesh % MeshDim
    NumberOfElements = Mesh % NumberOfBulkElements
    
    IF(.NOT. ASSOCIATED( Particles % InternalElements )) THEN
      ALLOCATE( Particles % InternalElements(NumberOfElements),STAT=istat )
      IF( istat /= 0 ) THEN
        CALL Fatal('MarkInternalElements','Allocation error 1')
      END IF
    END IF
    
    InternalElements => Particles % InternalElements
    InternalElements = .TRUE.
    
    DO t=1,NumberOfElements
      
      BulkElement => Mesh % Elements(t)
      
      body_id = BulkElement % BodyId
      
      IF(.FALSE.) THEN
        Body => CurrentModel % Bodies(body_id) % Values
        mat_id = ListGetInteger( Body,'Material',Found)
        bf_id = ListGetInteger( Body,'Body Force',Found)
      END IF
      
      IF( dim == 3 ) THEN
        imax = BulkElement % TYPE % NumberOfFaces 
      ELSE
        imax = BulkElement % TYPE % NumberOfEdges  
      END IF
      
      Hit = .FALSE.
      
      DO i=1, imax
        IF( dim == 3 ) THEN
          j = BulkElement % FaceIndexes(i)
          BoundaryElement => Mesh % Faces( j )
        ELSE
          j = BulkElement % EdgeIndexes(i)
          BoundaryElement => Mesh % Edges(j)
        END IF
        
        IF( .NOT. ASSOCIATED( BoundaryElement % BoundaryInfo ) ) CYCLE
        
        IF( ASSOCIATED( BulkElement, BoundaryElement % BoundaryInfo % Right ) ) THEN
          BulkElement2 => BoundaryElement % BoundaryInfo % Left
        ELSE
          BulkElement2 => BoundaryElement % BoundaryInfo % Right
        END IF
        
        ! A true boundary element
        IF( .NOT. ASSOCIATED( BulkElement2 )) THEN
          Hit = .TRUE.
          EXIT
        END IF
        
        body_id2 = BulkElement2 % BodyId
        IF( body_id2 == body_id ) CYCLE
        
        ! If the bodies are the same then there is no boundary
        IF(.TRUE.) THEN
          IF( body_id2 /= body_id ) THEN
            Hit = .TRUE.
            EXIT
          END IF
        ELSE        
          Body2 => CurrentModel % Bodies(body_id2) % Values
          
          mat_id2 = ListGetInteger( Body2,'Material')
          IF( mat_id2 /= mat_id ) THEN
            Hit = .TRUE.
            EXIT
          END IF
          
          bf_id2 = ListGetInteger( Body2,'Body Force',Found)
          IF( bf_id2 /= bf_id ) THEN
            Hit = .TRUE.
            EXIT
          END IF
        END IF
      END DO
      
      IF( Hit ) InternalElements(t) = .FALSE.
    END DO
    
    i = COUNT( InternalElements )
    j = NumberOfElements - i
    
    i = ParallelReduction( i ) 
    j = ParallelReduction( j ) 

    CALL Info('MarkInternalElements','Internal Elements: '//I2S(i),Level=8 )
    CALL Info('MarkInternalElements','Interface Elements: '//I2S(j),Level=8 )    
    
  END SUBROUTINE MarkInternalElements
  



  !---------------------------------------------------------
  !> Subroutine sets up some preliminary information needed for the 
  !> particler tracker: timeorder, space dimension, 
  !> bounding box, and mesh edges/faces.
  !---------------------------------------------------------
  SUBROUTINE SetParticlePreliminaries(Particles,dim,TimeOrder)
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: dim
    INTEGER, OPTIONAL :: TimeOrder
    
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: MinCoord(3), MaxCoord(3), s(3)
    INTEGER :: ierr
    
    Mesh => GetMesh()
    IF( .NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('SetParticleDimensions','No Mesh associated')
    END IF
    
    IF(PRESENT(TimeOrder)) THEN
      Particles % TimeOrder = TimeOrder
    ELSE
      Particles % TimeOrder = 2
    END IF
    
    IF( PRESENT( dim ) ) THEN
      IF( dim == 2 .OR. dim == 3 ) THEN
        Particles % dim = dim
      ELSE
        CALL Fatal('SetParticleDimensions','Invalid dimension')
      END IF
    ELSE
      Particles % dim = Mesh % Meshdim
    END IF
    
    MinCoord(1) = MINVAL(Mesh % Nodes % x )
    MinCoord(2) = MINVAL(Mesh % Nodes % y )
    MinCoord(3) = MINVAL(Mesh % Nodes % z )
    
    MaxCoord(1) = MAXVAL(Mesh % Nodes % x )
    MaxCoord(2) = MAXVAL(Mesh % Nodes % y )
    MaxCoord(3) = MAXVAL(Mesh % Nodes % z )
    
    Particles % LocalMinCoord = MinCoord
    Particles % LocalMaxCoord = MaxCoord
    
    
    ! Make a parallel reduction
    IF( ParEnv % PEs > 1 ) THEN
      s = MinCoord
      CALL MPI_ALLREDUCE( s, mincoord, 3, MPI_DOUBLE_PRECISION, &
          MPI_MIN, ELMER_COMM_WORLD, ierr )
      
      s = MaxCoord
      CALL MPI_ALLREDUCE( s, maxcoord, 3, MPI_DOUBLE_PRECISION, &
          MPI_MAX, ELMER_COMM_WORLD, ierr )
    END IF
    
    Particles % GlobalMinCoord = MinCoord
    Particles % GlobalMaxCoord = MaxCoord
        
    ! Create list of faces / edges 
    !-------------------------------------------------------------------------
    Mesh => GetMesh()
    CALL FindMeshEdges( Mesh )
    IF ( ParEnv % PEs > 1 ) THEN
      CALL SParEdgeNumbering(Mesh,Allmesh=.TRUE.)
      CALL SParFaceNumbering(Mesh,Allmesh=.TRUE.)
    END IF
    
    ! Mark elements that are not on boundary to make life faster in the future
    !-------------------------------------------------------------------------
    CALL MarkInternalElements( Particles )
    
  END SUBROUTINE SetParticlePreliminaries


  !---------------------------------------------------------
  !> Subroutine allocate particles before launching them.
  !---------------------------------------------------------
  SUBROUTINE AllocateParticles(Particles,NoParticles)
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: NoParticles
    
    REAL(KIND=dp), POINTER :: Velocity(:,:), Force(:,:), &
        Coordinate(:,:), PrevCoordinate(:,:), PrevVelocity(:,:)
    INTEGER, POINTER :: Status(:), ElementIndex(:), FaceIndex(:), NodeIndex(:), &
	Closest(:),Partition(:),Group(:)
    INTEGER :: PrevNoParticles, dofs, No, n, dim, TimeOrder, n1, n2, AllocParticles
    INTEGER, ALLOCATABLE :: Perm(:)
    
    IF( NoParticles <= Particles % MaxNumberOfParticles ) THEN
      CALL Info('AllocateParticles','There are already enough particles',Level=20)
      RETURN
    END IF

    PrevNoParticles = Particles % NumberOfParticles

    AllocParticles = NoParticles
    IF( PrevNoParticles > 0 ) THEN
      ! In parallel have a small buffer so that we are not next step here again!
      IF( ParEnv % PEs > 0 ) THEN
        AllocParticles = 1.02 * AllocParticles
      END IF
    END IF
          
    WRITE(Message,'(A,I0)') 'Allocating number of particles: ',AllocParticles
    CALL Info('AllocateParticles',Message,Level=12)    
    
    TimeOrder = Particles % TimeOrder
    dim = Particles % dim 
    dofs = dim

    ! Set pointers to the old stuff, these are needed
    ! if growing an already existing list of particles. 
    Coordinate => Particles % Coordinate
    PrevCoordinate => Particles % PrevCoordinate
    PrevVelocity => Particles % PrevVelocity
    Velocity => Particles % Velocity
    Force => Particles % Force
    Status => Particles % Status
    FaceIndex => Particles % FaceIndex
    ElementIndex => Particles % ElementIndex
    NodeIndex => Particles % NodeIndex
    Partition => Particles % Partition
    Group => Particles % Group

    ! Allocate the desired number of particles
    ALLOCATE( Particles % Coordinate(AllocParticles,dofs))
    ALLOCATE( Particles % Velocity(AllocParticles,dofs))
    ALLOCATE( Particles % Force(AllocParticles,dofs) )
    ALLOCATE( Particles % PrevCoordinate(AllocParticles,dofs) )
    
    ALLOCATE( Particles % Status(AllocParticles))
    ALLOCATE( Particles % ElementIndex(AllocParticles))
    ALLOCATE( Particles % FaceIndex(AllocParticles))

    IF( ASSOCIATED( PrevVelocity ) ) THEN
      ALLOCATE( Particles % PrevVelocity(AllocParticles,dofs) )
    END IF

    IF( Particles % NeighbourTable ) THEN
      Closest => Particles % ClosestNode        
      ALLOCATE( Particles % ClosestNode(AllocParticles) )
    END IF
    
    IF( Particles % NumberOfGroups > 0 ) THEN
      ALLOCATE( Particles % Group( AllocParticles ) )
    END IF
    
    IF( ASSOCIATED( NodeIndex ) ) THEN
      ALLOCATE( Particles % NodeIndex(AllocParticles) )
    END IF

    IF( ASSOCIATED( Partition ) ) THEN
      ALLOCATE( Particles % Partition(AllocParticles) )
      Particles % Partition = -1
    END IF
      
    ! Delete lost particles and move the remaining ones so that each 
    ! empty slot is used. Perm is the reordering of the particles, 
    ! not the reordering of field values as normally. 
    ! ---------------------------------------------------------------
    IF( PrevNoParticles > 0 ) THEN
      n = 0
      ALLOCATE( Perm( PrevNoParticles ) )
      Perm = 0
           
      DO No=1,PrevNoParticles
        IF ( Status(No) == PARTICLE_LOST ) CYCLE
        n = n+1
        Perm(n) = No
      END DO

      WRITE(Message,'(A,I0)') 'Number of old active particles: ',n
      CALL Info('AllocateParticles',Message,Level=8)    

      IF( n < PrevNoParticles ) THEN
        WRITE(Message,'(A,I0)') 'Number of deleted particles: ',PrevNoParticles-n
        CALL Info('AllocateParticles',Message,Level=10)    
      END IF

      n1 = 1
      n2 = n
      
      Particles % Coordinate(n1:n2,:) = Coordinate(Perm(n1:n2),:)
      Particles % Velocity(n1:n2,:) = Velocity(Perm(n1:n2),:)
      Particles % Force(n1:n2,:) = Force(Perm(n1:n2),:)
      Particles % PrevCoordinate(n1:n2,:) = PrevCoordinate(Perm(n1:n2),:)

      Particles % Status(n1:n2) = Status(Perm(n1:n2))
      Particles % FaceIndex(n1:n2) = FaceIndex(Perm(n1:n2))
      Particles % ElementIndex(n1:n2) = ElementIndex(Perm(n1:n2))

      IF( ASSOCIATED( PrevVelocity ) ) &
        Particles % PrevVelocity(n1:n2,:) = PrevVelocity(Perm(n1:n2),:)
      IF( Particles % NeighbourTable ) &
          Particles % ClosestNode(n1:n2) = Closest(Perm(n1:n2))
      IF ( ASSOCIATED(NodeIndex) ) &
          Particles % NodeIndex(n1:n2) = NodeIndex(Perm(n1:n2))
      IF ( ASSOCIATED(Partition) ) &
          Particles % Partition(n1:n2) = Partition(Perm(n1:n2))

      IF ( ASSOCIATED(Group) ) &
          Particles % Group(n1:n2) = Group(Perm(n1:n2))      

      PrevNoParticles = n
      Particles % NumberOfParticles = n
      
      ! Deallocate the old stuff
      DEALLOCATE( Coordinate, Velocity, Force, PrevCoordinate )
      DEALLOCATE( Status, FaceIndex, ElementIndex ) 

      IF( ASSOCIATED( PrevVelocity ) ) DEALLOCATE( PrevVelocity )
      IF( Particles % NeighbourTable ) DEALLOCATE(Closest)
      IF ( ASSOCIATED(NodeIndex) ) DEALLOCATE(NodeIndex)
      IF ( ASSOCIATED(Partition) ) DEALLOCATE(Partition)
      IF( ASSOCIATED( Group ) ) DEALLOCATE( Group ) 
    END IF

    ! Initialize the newly allocated particles with default values
    !-------------------------------------------------------------
    n1 = PrevNoParticles+1
    n2 = AllocParticles 

    Particles % Coordinate(n1:n2,:) = 0.0_dp
    Particles % Velocity(n1:n2,:) = 0.0_dp    
    Particles % Force(n1:n2,:) = 0.0_dp
    Particles % PrevCoordinate(n1:n2,:) = 0.0_dp

    Particles % Status(n1:n2) = PARTICLE_ALLOCATED
    Particles % ElementIndex(n1:n2) = 0
    Particles % FaceIndex(n1:n2) = 0
    
    IF( ASSOCIATED( Particles % PrevVelocity ) ) &
      Particles % PrevVelocity(n1:n2,:) = 0.0_dp

    IF( Particles % NeighbourTable ) &
        Particles % ClosestNode(n1:n2) = 0
    
    IF( ASSOCIATED( Particles % NodeIndex) ) &
        Particles % NodeIndex(n1:n2) = 0
    
    IF( ASSOCIATED( Particles % Partition) ) &
        Particles % Partition(n1:n2) = ParEnv % MyPe + 1

    IF( ASSOCIATED( Particles % Group ) ) &
        Particles % Group(n1:n2) = 0

    
    Particles % MaxNumberOfParticles = AllocParticles

    ! Finally resize the generic variables related to the particles
    !--------------------------------------------------------------
    CALL ParticleVariablesResize( Particles, PrevNoParticles, AllocParticles, Perm )

    IF( PrevNoParticles > 0 ) THEN
      CALL Info('AllocateParticles','Deallocating particle permutation',Level=20)
      DEALLOCATE( Perm )
    END IF

  END SUBROUTINE AllocateParticles

  !----------------------------------------------------
  !> Subroutine deletes lost particles. The reason for losing
  !> particles may be that they go to a neighbouring partition.
  !----------------------------------------------------
  SUBROUTINE DeleteLostParticles(Particles)
    TYPE(Particle_t), POINTER :: Particles
    
    INTEGER :: No, n, PrevNoParticles, n1, n2
    INTEGER, ALLOCATABLE :: Perm(:)

    PrevNoParticles = Particles % NumberOfParticles
    IF( PrevNoParticles == 0 ) RETURN

    ALLOCATE( Perm( PrevNoParticles ) )
    Perm = 0

    n = 0
    n1 = 0 
    DO No=1,PrevNoParticles
      IF ( Particles % Status(No) == PARTICLE_LOST ) CYCLE
      n = n + 1
      IF( n1 == 0 .AND. n /= No ) n1 = n
      Perm(n) = No
    END DO

    n2 = n
    CALL Info('DeleteLostParticles','Number of active particles: '&
        //I2S(n2),Level=12)

    IF(n1 == 0 ) THEN
      CALL Info('DeleteLostParticles','No particles need to be deleted',Level=12)
      RETURN
    ELSE
      CALL Info('DeleteLostParticles','First particle with changed permutation: '&
          //I2S(n1),Level=12)      
    END IF

    Particles % Coordinate(n1:n2,:) = &
        Particles % Coordinate(Perm(n1:n2),:)
    Particles % Velocity(n1:n2,:) = Particles % Velocity(Perm(n1:n2),:)
    Particles % Force(n1:n2,:) = Particles % Force(Perm(n1:n2),:)
    Particles % PrevCoordinate(n1:n2,:) = &
        Particles % PrevCoordinate(Perm(n1:n2),:)

    Particles % Status(n1:n2) = Particles % Status(Perm(n1:n2))
    Particles % FaceIndex(n1:n2) = Particles % FaceIndex(Perm(n1:n2))
    Particles % ElementIndex(n1:n2) = Particles % ElementIndex(Perm(n1:n2))

    IF( ASSOCIATED( Particles % PrevVelocity ) ) &
        Particles % PrevVelocity(n1:n2,:) = Particles % PrevVelocity(Perm(n1:n2),:)
    IF( Particles % NeighbourTable ) &
        Particles % ClosestNode(n1:n2) = Particles % ClosestNode(Perm(n1:n2))
    IF ( ASSOCIATED(Particles % NodeIndex) ) &
        Particles % NodeIndex(n1:n2) = Particles % NodeIndex(Perm(n1:n2))
    IF ( ASSOCIATED(Particles % Partition) ) &
        Particles % Partition(n1:n2) = Particles % Partition(Perm(n1:n2))
    
    IF( Particles % NumberOfGroups > 0 ) &
        Particles % Group(n1:n2) = Particles % Group(Perm(n1:n2))
    
    Particles % NumberOfParticles = n2
    
    IF ( n2 < PrevNoParticles ) THEN
       Particles % Coordinate(n2+1:PrevNoParticles,:) = 0._dp
       Particles % Velocity(n2+1:PrevNoParticles,:) = 0.0_dp
       Particles % Force(n2+1:PrevNoParticles,:) = 0._dp
       Particles % PrevCoordinate(n2+1:PrevNoParticles,:) = 0._dp      
      
       Particles % Status(n2+1:PrevNoParticles) = PARTICLE_ALLOCATED
       Particles % FaceIndex(n2+1:PrevNoParticles) = 0
       Particles % ElementIndex(n2+1:PrevNoParticles) = 0
       
       IF( ASSOCIATED( Particles % PrevVelocity ) ) &
           Particles % PrevVelocity(n2+1:PrevNoParticles,:) = 0._dp
       IF( Particles % NeighbourTable ) &
         Particles % ClosestNode(n2+1:PrevNoParticles) = 0
       IF ( ASSOCIATED(Particles % NodeIndex) ) &
          Particles % NodeIndex(n2+1:PrevNoParticles) = 0
       IF ( ASSOCIATED(Particles % Partition) ) &
           Particles % Partition(n2+1:PrevNoParticles) = 0      
       IF ( Particles % NumberOfGroups > 0 ) &
           Particles % Group(n2+1:PrevNoParticles) = 0      
    END IF

    ! Rorders the variables accordingly to Perm, no resize is needed since number of
    ! particles is not changed.
    CALL ParticleVariablesResize( Particles, PrevNoParticles, PrevNoParticles, Perm )


  END SUBROUTINE DeleteLostParticles
  


  !----------------------------------------------------
  !> Subroutine sets particles that are still sitting on boundary to 
  !> leave the premises and call for their deletion.
  !----------------------------------------------------
  SUBROUTINE EliminateExitingParticles( Particles )
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: NoParticles, No, DeletedParticles, dim, CumDeleted = 0, LimDeleted
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: BoundaryElement
    REAL(KIND=dp) :: Dist, Coord(3), Normal(3), SqrtElementMetric
    REAL(KIND=dp), POINTER :: Basis(:)
    TYPE(Nodes_t) :: BoundaryNodes
    INTEGER :: n
    LOGICAL :: Stat, Visited = .FALSE.
    
    SAVE Visited, CumDeleted, BoundaryNodes, Basis

RETURN
    
    NoParticles = Particles % NumberOfParticles
    dim = Particles % dim
    Mesh => GetMesh()
    
    IF(.NOT. Visited ) THEN
      Visited = .TRUE.
      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n) )
    END IF

    ! Just set the wall particles to lost elements for now    
    IF(.TRUE.) THEN
      DO No=1, NoParticles
        IF( Particles % Status(No) == PARTICLE_WALLBOUNDARY ) THEN
          Particles % Status(No) = PARTICLE_LOST          
        END IF
      END DO
    ELSE
    
      ! Some heuristics is used when to activate the deleting of particles
      ! Check for co-operation with parallel stuff.
      !-------------------------------------------------------------------
      LimDeleted = MAX( NINT( SQRT( 1.0_dp * NoParticles ) ), 20 )
      !  LimDeleted = 0
      
      DeletedParticles = 0
      DO No=1, NoParticles
        IF( Particles % Status(No) == PARTICLE_WALLBOUNDARY ) THEN
          DeletedParticles = DeletedParticles + 1
          Particles % Status(No) = PARTICLE_LOST
        END IF
      END DO
      
      CumDeleted = CumDeleted + DeletedParticles
      
      IF( CumDeleted > LimDeleted ) THEN
        !PRINT *,'Number of deleted particles:',CumDeleted, DeletedParticles, LimDeleted
        CALL DeleteLostParticles( Particles ) 
        CumDeleted = 0
      END IF
    END IF

  END SUBROUTINE EliminateExitingParticles
  

  !----------------------------------------------------
  !> Increase particle array size by given amount.
  !----------------------------------------------------
  SUBROUTINE IncreaseParticles(Particles,NoParticles)
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: NoParticles
    
    INTEGER :: Maxn
    
    ! Garbage collection:
    ! -------------------
    CALL DeleteLostParticles(Particles)
    
    ! Check if really need to allocate more space:
    ! --------------------------------------------
    Maxn = Particles % NumberOfParticles+NoParticles
    IF ( Maxn > Particles % MaxNumberOfParticles ) &
        CALL AllocateParticles( Particles, Maxn )
    
  END SUBROUTINE IncreaseParticles
  
  

  SUBROUTINE DestroyParticles(Particles) 
    TYPE(Particle_t), POINTER :: Particles

    CALL Info('DestroyParticles','Destrying the particle structures',Level=10)

    
    IF ( ASSOCIATED(Particles % Velocity) ) &
        DEALLOCATE( Particles % Velocity ) 
    
    IF ( ASSOCIATED(Particles % Force) ) &
        DEALLOCATE( Particles % Force )

    IF ( ASSOCIATED(Particles % PrevCoordinate) ) &
        DEALLOCATE( Particles % PrevCoordinate )

    IF ( ASSOCIATED(Particles % PrevVelocity) ) &
        DEALLOCATE( Particles % PrevVelocity )

    IF ( ASSOCIATED(Particles % NodeIndex) ) &
        DEALLOCATE( Particles % NodeIndex )

    IF ( ASSOCIATED(Particles % Partition) ) &
        DEALLOCATE( Particles % Partition )

    IF ( Particles % NumberOfGroups > 0 ) &
        DEALLOCATE( Particles % Group )

    IF( ASSOCIATED( Particles % Coordinate ) ) &
        DEALLOCATE( Particles % Coordinate ) 

    IF( ASSOCIATED( Particles % Status ) ) &
        DEALLOCATE( Particles % Status ) 

    IF( ASSOCIATED( Particles % FaceIndex ) ) &
        DEALLOCATE( Particles % FaceIndex )

    IF( ASSOCIATED( Particles % ElementIndex ) ) &
        DEALLOCATE( Particles % ElementIndex ) 
    
    IF( ASSOCIATED( Particles % UVW ) ) &
        DEALLOCATE( Particles % UVW )

    Particles % NumberOfParticles = 0
    Particles % MaxNumberOfParticles = 0
    
  END SUBROUTINE DestroyParticles
  


  !---------------------------------------------------------
  !> Subroutine for releaseing initiated but waiting particles.
  !---------------------------------------------------------
  SUBROUTINE ReleaseWaitingParticles(Particles) 
    TYPE(Particle_t), POINTER :: Particles
    
    TYPE(ValueList_t), POINTER :: Params
    INTEGER, POINTER :: Status(:)
    INTEGER :: i,j,NoParticles,ReleaseCount=0,ReleaseSet
    REAL(KIND=dp) :: ReleaseFraction
    LOGICAL :: Found,Visited = .FALSE.
    
    SAVE Visited, ReleaseCount
    
    ! Check whether all particles have already been released
    !-------------------------------------------------------
    NoParticles = Particles % NumberOfParticles
    IF( ReleaseCount >= NoParticles ) RETURN
    
    
    ! Get the size of the current release set
    !-------------------------------------------------------
    Params => ListGetSolverParams()
    ReleaseSet = GetInteger( Params,'Particle Release Number',Found)
    IF( .NOT. Found ) THEN
      ReleaseFraction = GetCReal( Params,'Particle Release Fraction',Found )
      IF(.NOT. Found ) THEN
        RETURN
      ELSE
        ReleaseSet = NINT( ReleaseFraction * NoParticles ) 
      END IF
    END IF
    CALL Info('ReleaseWaitingParticles','Releasing number of particles: '&
        //I2S(ReleaseCount),Level=10)
    
    IF( ReleaseSet <= 0 ) RETURN
    
    ! Release some waiting particles
    !-------------------------------------------------------
    Status => Particles % Status
    j = 0
    DO i=1,NoParticles
      IF( Status(i) == PARTICLE_WAITING ) THEN      
        Status(i) = PARTICLE_INITIATED 
        j = j + 1
        IF( j == ReleaseSet ) EXIT
      END IF
    END DO
    ReleaseCount = ReleaseCount + j
    
    
  END SUBROUTINE ReleaseWaitingParticles
  

  !---------------------------------------------------------
  !> Subroutine for chanching the partition of particles that
  !> cross the partition boundary.
  !---------------------------------------------------------
  FUNCTION ChangeParticlePartition(Particles) RESULT(nReceived)
    !---------------------------------------------------------
    TYPE(Particle_t), POINTER :: Particles
    !---------------------------------------------------------
    TYPE(Element_t), POINTER :: Face, Parent, Faces(:)
    
    INTEGER i,j,k,l,m,n,dim,NoPartitions, nextPart, nFaces, &
        Proc, ierr, status(MPI_STATUS_SIZE), n_part, nReceived, nSent, &
        ncomp, ncompInt
    
    INTEGER, ALLOCATABLE :: Perm(:), Indexes(:), Neigh(:), &
        Recv_parts(:), Requests(:)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var

    LOGICAL :: Found
    INTEGER, POINTER :: Neighbours(:)
    LOGICAL, POINTER :: FaceInterface(:), IsNeighbour(:)

    INTEGER :: q, CntFailed
    LOGICAL, ALLOCATABLE :: Failed(:)
    TYPE(ValueList_t), POINTER :: BC
    INTEGER, ALLOCATABLE :: BCcount(:)
    
    TYPE ExchgInfo_t
      INTEGER :: n=0
      INTEGER, ALLOCATABLE :: Gindex(:), Particles(:)
    END TYPE ExchgInfo_t
    
    REAL(KIND=dp), ALLOCATABLE :: Buf(:)
    INTEGER, ALLOCATABLE :: BufInt(:)
    TYPE(ExchgInfo_t), ALLOCATABLE :: ExcInfo(:)
    CHARACTER(*), PARAMETER :: Caller = 'ChangeParticlePartition'
    !---------------------------------------------------------
    
    nReceived = 0
    IF( ParEnv% PEs == 1 ) RETURN
    
    CALL Info(Caller,'Sending particles among partitions',Level=10)

    Mesh => GetMesh()
    dim = Particles % dim
    
    ! Count & Identify neighbouring partitions:
    ! -----------------------------------------
    ALLOCATE(IsNeighbour(ParEnv % PEs))
    NoPartitions = MeshNeighbours(Mesh,IsNeighbour)
    ALLOCATE(Perm(ParEnv % PEs), Neigh(NoPartitions) )
    Perm = 0
    
    ALLOCATE( BCCount(CurrentModel % NumberOfBCs) )
    BCCount = 0

    CntFailed = 0
    
    NoPartitions=0
    DO i=1,ParEnv % PEs
      IF ( i-1 == ParEnv % Mype ) CYCLE
      IF ( IsNeighbour(i) ) THEN
        NoPartitions = NoPartitions+1
        Perm(i) = NoPartitions
        Neigh(NoPartitions) = i-1
      END IF
    END DO
    DEALLOCATE(IsNeighbour)

    CALL Info(Caller,'Number of neighbour partitions: '&
        //I2S(NoPartitions),Level=12)
    
    !
    ! Count particles to be sent to neighbours:
    ! -----------------------------------------
    ALLOCATE(ExcInfo(NoPartitions))
    ExcInfo % n = 0
    
    DO i=1,Particles % NumberOfParticles
      IF( Particles % Status(i) /= PARTICLE_WALLBOUNDARY ) CYCLE
      
      IF ( dim==2 ) THEN
        Face => Mesh % Edges(Particles % FaceIndex(i))
        FaceInterface  => Mesh % ParallelInfo % EdgeInterface
        Neighbours => Mesh % ParallelInfo %  &
            EdgeNeighbourList(Face % ElementIndex) % Neighbours
      ELSE
        Face => Mesh % Faces(Particles % FaceIndex(i))
        FaceInterface => Mesh % ParallelInfo % FaceInterface
        Neighbours => Mesh % ParallelInfo %  &
            FaceNeighbourList(Face % ElementIndex) % Neighbours
      END IF
      
      IF ( FaceInterface(Face % ElementIndex) ) THEN
        BC => NULL()
        IF( ASSOCIATED( Face % BoundaryInfo ) ) THEN
          k = Face % BoundaryInfo % Constraint
          IF( k > 0 ) THEN
            DO j=1,CurrentModel % NumberOfBCs
              IF ( k == CurrentModel % BCs(j) % Tag ) THEN
                BC => CurrentModel % BCs(j) % Values
                EXIT
              END IF
            END DO
          END IF
        END IF
        IF( ASSOCIATED( BC ) ) THEN
          BCCount(j) = BCCount(j) + 1          
          IF( ListGetLogical( BC,'Particle Wall',Found) ) THEN
            Particles % Status(i) = PARTICLE_WALLBOUNDARY
            CYCLE
          END IF
        END IF
        
        nextPart = ParEnv % MyPE
        DO j=1,SIZE(Neighbours)
          IF ( ParEnv % Mype /= Neighbours(j) ) THEN
            nextPart = Neighbours(j)
            k = Perm(nextPart+1)
            IF ( k>0 ) THEN
              ExcInfo(k) % n = ExcInfo(k) % n+1
              Particles % Status(i) = PARTICLE_PARTBOUNDARY
            ELSE 
              Particles % Status(i) = PARTICLE_LOST
            END IF
            EXIT
          END IF
        END DO
      END IF
    END DO

    DO j=1,CurrentModel % NumberOfBCs
      IF( BCCount(j) > 0 ) THEN
        CALL Warn(Caller,I2S(BCCount(j))//' particles hit BC '//I2S(j))
      END IF
    END DO

    
    n = SUM( ExcInfo(1:NoPartitions) % n )
    CALL Info(Caller,'Number of particles to send: '&
        //I2S(n),Level=10)
    
    CALL MPI_ALLREDUCE( n, nSent, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )
    IF ( nSent == 0 ) THEN
      CALL Info(Caller,'No particles needs to be sent',Level=10)
      DEALLOCATE(ExcInfo, Perm, Neigh)
      RETURN
    ELSE
      CALL Info(Caller,'Global number of particles to sent: '&
          //I2S(nSent),Level=10)      
    END IF

    !
    ! Receive interface sizes:
    !--------------------------
    ALLOCATE( Recv_Parts(NoPartitions), Requests(NoPartitions) )
    DO i=1,NoPartitions
      CALL MPI_iRECV( Recv_Parts(i),1, MPI_INTEGER, Neigh(i), &
          1000, ELMER_COMM_WORLD, requests(i), ierr )
    END DO
    
    DO i=1,NoPartitions
      CALL MPI_BSEND( ExcInfo(i) % n, 1, MPI_INTEGER, Neigh(i), &
          1000, ELMER_COMM_WORLD, ierr )
    END DO
    CALL MPI_WaitAll( NoPartitions, Requests, MPI_STATUSES_IGNORE, ierr )

    n = SUM(Recv_Parts)
    CALL Info(Caller,'Number of particles to receive: '&
        //I2S(n),Level=10)

    CALL MPI_ALLREDUCE( n, nReceived, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )
    IF ( nReceived==0 ) THEN
      CALL Info(Caller,'No particles needs to be received',Level=10)
      DEALLOCATE(Recv_Parts, Requests, ExcInfo, Perm, Neigh)
      RETURN
    ELSE
      CALL Info(Caller,'Global number of particles to receive: '&
          //I2S(nReceived),Level=10)
    END IF

    n = SUM( ExcInfo(1:NoPartitions) % n )
    CALL Info(Caller,'Total number of particles to sent: '&
        //I2S(n),Level=10)


    !
    ! Collect particles to be sent to neighbours:
    ! -------------------------------------------
    DO i=1,NoPartitions
      ALLOCATE( ExcInfo(i) % Gindex(ExcInfo(i) % n), &
          ExcInfo(i) % Particles(ExcInfo(i) % n) )
      ExcInfo(i) % n = 0
    END DO
    
    DO i=1,Particles % NumberOfParticles
      IF( Particles % Status(i) /= PARTICLE_PARTBOUNDARY ) CYCLE
      
      IF ( dim==2 ) THEN
        Face => Mesh % Edges(Particles % FaceIndex(i))
        FaceInterface  => Mesh % ParallelInfo % EdgeInterface
        Neighbours => Mesh % ParallelInfo %  &
            EdgeNeighbourList(Face % ElementIndex) % Neighbours
      ELSE
        Face => Mesh % Faces(Particles % FaceIndex(i))
        FaceInterface => Mesh % ParallelInfo % FaceInterface
        Neighbours => Mesh % ParallelInfo %  &
            FaceNeighbourList(Face % ElementIndex) % Neighbours
      END IF
      
      IF ( FaceInterface(Face % ElementIndex) ) THEN
        nextPart = ParEnv % MyPE
        DO j=1,SIZE(Neighbours)
          IF ( ParEnv % Mype /= Neighbours(j) ) THEN
            nextPart = Neighbours(j);
            EXIT
          END IF
        END DO
        Particles % Status(i) = PARTICLE_LOST
        j = Perm(nextPart+1)
        IF ( j==0 ) THEN
          CALL Warn( Caller, 'Neighbouring partition not found?')
          CYCLE
        END IF
        ExcInfo(j) % n = ExcInfo(j) % n+1
        n = ExcInfo(j) % n
        ExcInfo(j) % Particles(n) = i
        ExcInfo(j) % Gindex(n) = Face % GElementIndex
      END IF
    END DO
    
    n = 0
    DO i=1,NoPartitions
      n = n + ExcInfo(i) % n
    END DO

    CALL Info(Caller,'Collected particles from partitions: '&
        //I2S(n),Level=12)

    ncomp = dim ! coordinate
    IF( ASSOCIATED( Particles % Velocity ) ) ncomp = ncomp + dim
    IF( ASSOCIATED( Particles % PrevCoordinate ) ) ncomp = ncomp + dim
    IF( ASSOCIATED( Particles % PrevVelocity ) ) ncomp = ncomp + dim
    IF( ASSOCIATED( Particles % Force ) ) ncomp = ncomp + dim
    Var => Particles % Variables
    DO WHILE( ASSOCIATED(Var) )
      IF( Var % Dofs == 1 ) THEN
        ncomp = ncomp + Var % Dofs 
      ELSE
        CALL Warn(Caller,'Implement for vectors only!')
      END IF
      Var => Var % Next 
    END DO

    ncompInt = 0
    IF ( ASSOCIATED(Particles % NodeIndex) ) ncompInt = ncompInt + 1
    IF ( ASSOCIATED(Particles % Partition) ) ncompInt = ncompInt + 1
    ! status, elementindex & closestnode are recomputed, and hence not communicated

    CALL Info(Caller,'Transferring real entries between particles: '&
        //I2S(ncomp),Level=12)
    CALL Info(Caller,'Transferring integer entries between particles: '&
        //I2S(ncompInt),Level=12)


    n = 2*(n + 2*ncomp + MPI_BSEND_OVERHEAD*2*NoPartitions)
    CALL CheckBuffer(n)
    
    CALL Info(Caller,'Size of data buffer: ' &
        //I2S(n),Level=12)

    ! Send particles:
    ! ---------------
    CALL Info(Caller,'Now sending particle data',Level=14)
    DO j=1,NoPartitions
      n = ExcInfo(j) % n
      IF ( n<=0 ) CYCLE
      
      CALL MPI_BSEND( ExcInfo(j) % Gindex, n, MPI_INTEGER, Neigh(j), &
          1001, ELMER_COMM_WORLD, ierr )
      
      ALLOCATE(Buf(dim*n*ncomp))
      m = 0
      DO k=1,dim
        DO l=1,n
          m = m + 1
          Buf(m) = Particles % Coordinate(ExcInfo(j) % Particles(l),k)
        END DO
      END DO
      
      IF ( ASSOCIATED(Particles % Velocity) ) THEN
        DO k=1,dim
          DO l=1,n
            m = m + 1
            Buf(m) = Particles % Velocity(ExcInfo(j) % Particles(l),k)
          END DO
        END DO
      END IF
      
      IF ( ASSOCIATED(Particles % Force) ) THEN
        DO k=1,dim
          DO l=1,n
            m = m + 1
            Buf(m) = Particles % Force(ExcInfo(j) % Particles(l),k)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED(Particles % PrevCoordinate) ) THEN
        DO k=1,dim
          DO l=1,n
            m = m + 1
            Buf(m) = Particles % PrevCoordinate(ExcInfo(j) % Particles(l),k)
          END DO
        END DO
      END IF
      
      IF ( ASSOCIATED(Particles % PrevVelocity) ) THEN
        DO k=1,dim
          DO l=1,n
            m = m + 1
            Buf(m) = Particles % PrevVelocity(ExcInfo(j) % Particles(l),k)
          END DO
        END DO
      END IF
      
      Var => Particles % Variables
      DO WHILE( ASSOCIATED(Var) )
        IF( Var % Dofs == 1 ) THEN
          DO l=1,n
            m = m + 1
            Buf(m) = Var % Values(ExcInfo(j) % Particles(l))
          END DO
        END IF
        Var => Var % Next 
      END DO
     
      CALL MPI_BSEND( Buf, m, MPI_DOUBLE_PRECISION, &
          Neigh(j), 1002, ELMER_COMM_WORLD, ierr )
      
      DEALLOCATE(Buf)

      IF( ncompInt >  0 ) THEN
        ! Sent integers also 
        ALLOCATE(BufInt(n*ncompInt))
        
        m = 0
        IF ( ASSOCIATED(Particles % NodeIndex) ) THEN
          DO l=1,n
            m = m + 1
            BufInt(m) = Particles % NodeIndex(ExcInfo(j) % Particles(l))
          END DO
        END IF
        
        IF ( ASSOCIATED(Particles % Partition) ) THEN
          DO l=1,n
            m = m + 1
            BufInt(m) = Particles % Partition(ExcInfo(j) % Particles(l))
         END DO
        END IF
        CALL MPI_BSEND( BufInt, m, MPI_INTEGER, Neigh(j), 1003, ELMER_COMM_WORLD, ierr )
        
        DEALLOCATE(BufInt)
      END IF
    END DO
    

    DEALLOCATE(Perm)
    DO i=1,NoPartitions
      DEALLOCATE( ExcInfo(i) % Gindex, ExcInfo(i) % Particles )
    END DO
    DEALLOCATE(ExcInfo)

    ! Recv particles:
    ! ---------------
    CALL Info(Caller,'Now receiving particle data',Level=14)

    n = SUM(Recv_Parts)

    CALL DeleteLostParticles(Particles)
    IF ( Particles % NumberOfParticles+n > Particles % MaxNumberOfParticles ) THEN
      CALL IncreaseParticles( Particles, Particles % NumberOfParticles + 2*n - &
                    Particles % MaxNumberOfParticles )
    END IF
    
    IF(Particles % dim==2 ) THEN
      nFaces = Mesh % NumberOfEdges
      Faces => Mesh % Edges
    ELSE
      Faces => Mesh % Faces
      nFaces = Mesh % NumberOfFaces
    END IF

    DO i=1,NoPartitions
      n = Recv_Parts(i)
      IF ( n<=0 ) CYCLE

      CALL Info(Caller,'Partition '//I2S(ParEnv % MyPe+1)//&
          ' receiving from '//I2S(i)//': '//TRIM(I2S(n)),Level=20)   
      
      proc = Neigh(i)
      ALLOCATE(Indexes(n))
      
      CALL MPI_RECV( Indexes, n, MPI_INTEGER, proc, &
          1001, ELMER_COMM_WORLD, status, ierr )

      ALLOCATE(Failed(n))
      Failed = .FALSE.
      
      n_part=Particles % NumberOfParticles
      DO j=1,n
        k=SearchElement( nFaces, Faces, Indexes(j) )
        IF ( k <= 0 ) THEN
          Failed(j) = .TRUE.
          CYCLE
        END IF
        
        Face => Faces(k) 
        Parent => Face % BoundaryInfo % Left
        IF ( .NOT.ASSOCIATED(Parent) ) &
            Parent => Face % BoundaryInfo % Right
        
        n_part = n_part+1
        Particles % Status(n_part) = PARTICLE_PARTBOUNDARY
        Particles % ElementIndex(n_part) = Parent % ElementIndex
      END DO
      
      m = n*ncomp*dim
      IF ( ASSOCIATED(Particles % Velocity) ) m=m+n*dim
      ALLOCATE(Buf(m))
      CALL MPI_RECV( Buf, m, MPI_DOUBLE_PRECISION, proc, &
          1002, ELMER_COMM_WORLD, status, ierr )
      
      n_part=Particles % NumberOfParticles
      m = 0
      DO k=1,dim
        q = 0
        DO l=1,n
          m = m + 1
          IF(Failed(l)) CYCLE
          q = q + 1
          Particles % Coordinate(n_part+q,k) = Buf(m)
        END DO
      END DO
      
      IF ( ASSOCIATED(Particles % Velocity) ) THEN
        DO k=1,dim
          q = 0
          DO l=1,n
            m = m + 1
            IF(Failed(l)) CYCLE
            q = q + 1
            Particles % Velocity(n_part+q,k) = Buf(m)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED(Particles % Force) ) THEN
        DO k=1,dim
          q = 0
          DO l=1,n
            m = m + 1
            IF(Failed(l)) CYCLE
            q = q + 1
            Particles % Force(n_part+q,k) = Buf(m)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED(Particles % PrevCoordinate) ) THEN
        DO k=1,dim
          q = 0
          DO l=1,n
            m = m + 1
            IF (Failed(l)) CYCLE
            q = q + 1
            Particles % PrevCoordinate(n_part+q,k) = Buf(m)
          END DO
        END DO
      END IF

      IF ( ASSOCIATED(Particles % PrevVelocity) ) THEN
        DO k=1,dim
          q  = 0
          DO l=1,n
            m = m + 1
            IF (Failed(l)) CYCLE
            q = q + 1
            Particles % PrevVelocity(n_part+q,k) = Buf(m)
          END DO
        END DO
      END IF

      Var => Particles % Variables
      DO WHILE( ASSOCIATED(Var) )
        IF( Var % Dofs == 1 ) THEN
          q = 0
          DO l=1,n
            m = m + 1
            IF(Failed(l)) CYCLE
            q = q + 1
            Var % Values(n_part+q) = Buf(m)
          END DO
        END IF
        Var => Var % Next 
      END DO

      DEALLOCATE(Buf)

      IF( ncompInt > 0 ) THEN

        ALLOCATE(BufInt(n*ncompInt))
        m = n*ncompInt
        
       CALL MPI_RECV( BufInt, m, MPI_INTEGER, proc, 1003, ELMER_COMM_WORLD, status, ierr )
        
       m = 0
       IF( ASSOCIATED( Particles % NodeIndex ) ) THEN
         q = 0
         DO l=1,n
           m = m + 1
           IF(Failed(l)) CYCLE
           q = q + 1
           Particles % NodeIndex(n_part+q) = BufInt(m)
         END DO
       END IF

       IF( ASSOCIATED( Particles % Partition ) ) THEN
         q = 0
         DO l=1,n
           m = m + 1
           IF(Failed(l)) CYCLE
           q = q + 1
           Particles % Partition(n_part+q) = BufInt(m)
         END DO
       END IF

       DEALLOCATE(BufInt)
      END IF

      CntFailed = CntFailed + COUNT(Failed)
      
      Particles % NumberOfParticles = Particles % NumberOfParticles + COUNT(.NOT.Failed)
      DEALLOCATE(Indexes, Failed)
    END DO

    IF( CntFailed > 0 ) THEN
      CALL Warn(Caller,'Particles not found on the other side: '//I2S(CntFailed))
    END IF
    
    DEALLOCATE(Recv_Parts, Neigh, Requests)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    CALL Info(Caller,'Information exchange done',Level=10)

    
  CONTAINS

   
    !
    ! Search an element Item from an ordered Element_t array(N) and return
    ! Index to that array element. Return value -1 means Item was not found.
    !
    FUNCTION SearchElement( N, IArray, Item ) RESULT(Indx)
      IMPLICIT NONE
      
      INTEGER :: Item, Indx, i
      INTEGER :: N
      TYPE(Element_t) :: Iarray(:)
      
      ! Local variables
      
      INTEGER :: Lower, Upper, lou
      
      !*********************************************************************
      
      Indx  = -1
      Upper =  N
      Lower =  1
      
      ! Handle the special case
      
      IF ( Upper < Lower ) RETURN
      
      DO WHILE( .TRUE. )
        IF ( IArray(Lower) % GelementIndex == Item ) THEN
          Indx = Lower
          EXIT
        ELSE IF ( IArray(Upper) % GelementIndex == Item ) THEN
          Indx = Upper
          EXIT
        END IF
        
        IF ( (Upper - Lower) > 1 ) THEN
          Lou = ISHFT((Upper + Lower), -1)
          IF ( IArray(lou) % GelementIndex < Item ) THEN
            Lower = Lou
          ELSE
            Upper = Lou
          END IF
        ELSE
          EXIT
        END IF
      END DO
    END FUNCTION SearchElement
    
  END FUNCTION ChangeParticlePartition
  

  !-----------------------------------------------------------------
  !> Subroutine for advecting back to given partition and node index.
  !-----------------------------------------------------------------
  SUBROUTINE ParticleAdvectParallel(Particles, SentField, RecvField, dofs ) 
    !----------------------------------------------------------------------
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp), POINTER :: SentField(:), RecvField(:)
    INTEGER :: Dofs
    !---------------------------------------------------------
    INTEGER i,j,k,idof,l,m,n,dim,NoPartitions,NoParticles,part, &
        ierr, status(MPI_STATUS_SIZE), nReceived, nSent, nerr
    INTEGER, ALLOCATABLE :: RecvParts(:), SentParts(:), Requests(:)
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), ALLOCATABLE :: SentReal(:),RecvReal(:)
    INTEGER, ALLOCATABLE :: SentInt(:),RecvInt(:)
    !---------------------------------------------------------
    
    CALL Info('ParticleAdvectParallel',&
        'Returning particle info to their initiating partition',Level=15)

    nReceived = 0
    IF( ParEnv% PEs == 1 ) RETURN
    
    Mesh => GetMesh()
    dim = Particles % dim

!debug = particles % partition(no) == 2 .and. particles % nodeindex(no) == 325
    NoPartitions = ParEnv % PEs
    NoParticles = Particles % NumberOfParticles

    IF( .NOT. ASSOCIATED( Particles % Partition ) ) THEN
      CALL Fatal('ParticleAdvectParallel','Partition must be present!')
    END IF
    IF( .NOT. ASSOCIATED( Particles % NodeIndex ) ) THEN
      CALL Fatal('ParticleAdvectParallel','NodeIndex must be present!')
    END IF

    ! First take the components of the field that lie in the present partition. 
    !-------------------------------------------------------------------------
    DO i=1,Particles % NumberOfParticles
      Part = Particles % Partition(i)
      IF( Part-1 == ParEnv % MyPe ) THEN
        j = Particles % NodeIndex(i)
        IF ( j==0 ) CYCLE
        IF( dofs == 1 ) THEN
          RecvField(j) = SentField(i)
        ELSE
          DO idof=1,dofs
            RecvField(dofs*(j-1)+idof) = SentField(dofs*(i-1)+idof)
          END DO
        END IF          
      END IF
    END DO

    ! Count particles to be sent to neighbours:
    ! -----------------------------------------
    ALLOCATE(SentParts(NoPartitions), RecvParts(NoPartitions), &
        Requests(NoPartitions))
    SentParts = 0
    RecvParts = 0
    Requests = 0

    nerr = 0
    DO i=1,Particles % NumberOfParticles
      Part = Particles % Partition(i)
      IF( Part-1 == ParEnv % MyPe ) CYCLE
      IF( Part < 1 .OR. Part > ParEnv % PEs ) THEN
        nerr = nerr + 1
        CYCLE
      END IF
      SentParts(Part) = SentParts(Part) + 1
    END DO
    IF( nerr > 0 ) THEN
      CALL Info('ParticleAdvectParallel','Invalid partition in particles: '//I2S(nerr))
    END IF
    
    n = SUM( SentParts )
    CALL Info('ParticleAdvectParallel','Local particles to be sent: '//I2S(n),Level=12)
    
    CALL MPI_ALLREDUCE( n, nSent, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )

    IF( nSent > 0 ) THEN
      CALL Info('ParticleAdvectParallel','Global particles to be sent: '&
          //I2S(nSent),Level=12)
    ELSE
      ! If nobody is sending any particles then there can be no need to receive particles either
      ! Thus we can make an early exit.
      DEALLOCATE(SentParts, RecvParts, Requests )
      CALL Info('ParticleAdvectParallel','Nothing to do in parallel!',Level=15)
      RETURN
    END IF

    ! Receive interface sizes:
    !--------------------------
    DO i=1,NoPartitions
      IF( i-1 == ParEnv % MyPe ) CYCLE
      CALL MPI_BSEND( SentParts(i), 1, MPI_INTEGER, i-1, &
          1000, ELMER_COMM_WORLD, ierr )
    END DO

    DO i=1,NoPartitions
      IF( i-1 == ParEnv % MyPe ) CYCLE
      CALL MPI_RECV( RecvParts(i), 1, MPI_INTEGER, i-1, &
          1000, ELMER_COMM_WORLD, Status, ierr )
    END DO

    n = SUM(RecvParts)
    CALL Info('ParticleAdvectParallel','Particles to be received: '//I2S(n),Level=12)

    CALL MPI_ALLREDUCE( n, nReceived, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )

    IF ( nReceived==0 ) THEN
      DEALLOCATE(RecvParts, SentParts, Requests )
      RETURN
    END IF

    CALL Info('ParticleAdvectParallel','Total number of particles to be received: '&
        //I2S(nReceived),Level=12)

    n = 2*(dofs*n + MPI_BSEND_OVERHEAD*2*NoPartitions)
    CALL CheckBuffer(n)

    CALL Info('ParticleAdvectParallel','Buffer size for sending and receiving: ' &
        //I2S(n),Level=14)


    ! Allocate sent and receive buffers based on the maximum needed size.
    !--------------------------------------------------------------------
    n = MAXVAL( SentParts )
    ALLOCATE( SentReal(dofs*n), SentInt(n) )
    CALL Info('ParticleAdvectParallel','Allocating sent buffer for size: '&
        //I2S(n),Level=18)

    n = MAXVAL( RecvParts )
    ALLOCATE( RecvReal(dofs*n), RecvInt(n) )
    CALL Info('ParticleAdvectParallel','Allocating receive buffer for size: '&
        //I2S(n),Level=18)

    ! Send particles:
    ! ---------------
    CALL Info('ParticleAdvectParallel','Now sending field',Level=14)

    DO j=1,NoPartitions

      IF( j-1 == ParEnv % MyPe ) CYCLE
      IF( SentParts(j) == 0 ) CYCLE

      m = 0
      DO l=1,NoParticles
        IF( Particles % Partition(l) == j ) THEN
          m = m + 1
          IF(dofs==1) THEN
            SentReal(m) = SentField(l)
          ELSE
            DO idof=1,dofs
              SentReal(dofs*(m-1)+idof) = SentField(dofs*(l-1)+idof)
            END DO
          END IF
          SentInt(m) = Particles % NodeIndex(l)
        END IF
      END DO

      CALL MPI_BSEND( SentInt, m, MPI_INTEGER, j-1, &
          1001, ELMER_COMM_WORLD, ierr )

      CALL MPI_BSEND( SentReal, dofs*m, MPI_DOUBLE_PRECISION, j-1, &
          1002, ELMER_COMM_WORLD, ierr )
    END DO


    ! Recv particles:
    ! ---------------
    CALL Info('ParticleAdvectParallel','Now receiving field',Level=14)

    nerr = 0
    DO j=1,NoPartitions

      IF( j-1 == ParEnv % MyPe ) CYCLE

      m = RecvParts(j)
      IF ( m == 0 ) CYCLE
      
      CALL MPI_RECV( RecvInt, m, MPI_INTEGER, j-1, &
          1001, ELMER_COMM_WORLD, status, ierr )

      CALL MPI_RECV( RecvReal, dofs*m, MPI_DOUBLE_PRECISION, j-1, &
          1002, ELMER_COMM_WORLD, status, ierr )
      
      DO l=1,m
        k = RecvInt(l)
        IF( k <=0 .OR. k > SIZE( RecvField ) ) THEN
          nerr = nerr + 1
          CYCLE
        END IF
        IF(dofs==1) THEN
          RecvField(k) = RecvReal(l)
        ELSE
          DO idof=1,dofs
            RecvField(dofs*(k-1)+idof) = RecvReal(dofs*(l-1)+idof)
          END DO
        END IF
      END DO
    END DO

    IF( nerr > 0 ) THEN
      CALL Info('ParticleAdvectParallel','Invalid received index in particles: '//I2S(nerr))
    END IF

    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    DEALLOCATE(SentParts, RecvParts, Requests, SentInt, SentReal, RecvInt, RecvReal )

    CALL Info('ParticleAdvectParallel','Particle field communication done',Level=14)

  END SUBROUTINE ParticleAdvectParallel




  !---------------------------------------------------------
  !> Subroutine computes the means of coordinates / velocities / force.
  ! The statistics could be made more detailed...
  !---------------------------------------------------------
  SUBROUTINE ParticleStatistics( Particles, DerOrder ) 
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: DerOrder
    
    REAL(KIND=dp) :: Coord(3),MeanCoord(3),AbsCoord(3),VarCoord(3), &
        MinCoord(3), MaxCoord(3),val    
    INTEGER :: i,j,Cnt,NoParticles,TotParticles,dim
    REAL(KIND=dp), POINTER :: TargetVector(:,:), StdVector(:)
    INTEGER, POINTER :: Status(:)
    CHARACTER(:), ALLOCATABLE :: DataName
    TYPE(Variable_t), POINTER :: PartVar
        
    MeanCoord = 0.0_dp
    AbsCoord = 0.0_dp
    VarCoord = 0.0_dp
    MinCoord = HUGE( MinCoord )
    MaxCoord = -HUGE( MaxCoord )
    
    Cnt = 0
    NoParticles =  Particles % NumberOfParticles
    dim = Particles % dim
    Coord = 0.0_dp
    StdVector => NULL()

    SELECT CASE( DerOrder )
    CASE(-1)
      TargetVector => Particles % PrevCoordinate
      DataName = 'previous coordinate values'
    CASE(0)
      TargetVector => Particles % Coordinate
      DataName = 'current coordinate values'      
    CASE(1)
      TargetVector => Particles % Velocity
      DataName = 'current velocity values'
    CASE(2)
      TargetVector => Particles % Force
      DataName = 'current force values'
    CASE(3)
      DataName = 'particle time integral'
      PartVar => ParticleVariableGet( Particles,DataName)
      IF(.NOT. ASSOCIATED(PartVar) ) RETURN
      StdVector => PartVar % Values 
    CASE(4)
      DataName = 'particle distance integral'
      PartVar => ParticleVariableGet( Particles,DataName)
      IF(.NOT. ASSOCIATED(PartVar) ) RETURN
      StdVector => PartVar % Values 
    CASE DEFAULT
      CALL Fatal('ParticleStatistics','Unknown value for DerOrder!')
    END SELECT
    
    Status => Particles % Status
    IF(ASSOCIATED(StdVector)) dim = 1
    
    DO i=1,NoParticles
      IF( Status(i) >= PARTICLE_LOST ) CYCLE
      IF( Status(i) < PARTICLE_INITIATED ) CYCLE

      IF( ASSOCIATED( StdVector ) ) THEN
        Coord(1) = StdVector(i)
      ELSE
        Coord(1:dim) = TargetVector(i,1:dim)
      END IF
        
      MeanCoord = MeanCoord + Coord
      AbsCoord = AbsCoord + ABS( Coord )
      VarCoord = VarCoord + Coord**2
      DO j=1,dim
        MinCoord(j) = MIN( MinCoord(j), Coord(j) )
        MaxCoord(j) = MAX( MaxCoord(j), Coord(j) )
      END DO
      Cnt = Cnt + 1
    END DO
    
    TotParticles = ParallelReduction( Cnt ) 
    IF( TotParticles == 0 ) THEN
      CALL Warn('MeanParticleCoordinate','No active particles!')
      RETURN
    END IF

    
    IF( TotParticles > Cnt ) THEN
      ! Compute parallel sums
      DO j=1,dim
        MeanCoord(j) = ParallelReduction( MeanCoord(j) )
        AbsCoord(j) = ParallelReduction( AbsCoord(j) ) 
        VarCoord(j) = ParallelReduction( varCoord(j) )
        MinCoord(j) = ParallelReduction( MinCoord(j),1 )
        MaxCoord(j) = ParallelReduction( MaxCoord(j),2 )
      END DO
    END IF

    IF( TotParticles == 1 ) THEN
      IF( ParEnv % myPE == 0 ) THEN       
        PRINT *,'Particle info on '//TRIM(DataName) 
        PRINT *,'Val: ',MinCoord(1:dim),&
            ' Abs: ',SQRT(SUM( MinCoord(1:dim)**2 ))
      END IF
    ELSE 
      MeanCoord = MeanCoord / TotParticles
      AbsCoord = AbsCoord / TotParticles
      
      ! If unlucky with rounding error this could be negative even though it really can't
      DO j=1,dim
        val = VarCoord(j) / TotParticles - MeanCoord(j)**2  
        IF( val <= TINY(val) ) THEN
          VarCoord(j) = 0.0_dp
        ELSE
          VarCoord(j) = SQRT( val )
        END IF
      END DO

      IF( ParEnv % myPE == 0 ) THEN       
        PRINT *,'Statistical info on '//TRIM(DataName) 
        PRINT *,'Mean:',MeanCoord(1:dim)
        PRINT *,'Abs: ',AbsCoord(1:dim)
        PRINT *,'Var: ',VarCoord(1:dim)
        PRINT *,'Min: ',MinCoord(1:dim)
        PRINT *,'Max: ',MaxCoord(1:dim)
      END IF
    END IF
    
  END SUBROUTINE ParticleStatistics
  

  !---------------------------------------------------------
  !> Echos some information on particle amounts to standard output.
  !---------------------------------------------------------
  SUBROUTINE ParticleInformation( Particles, ParticleStepsTaken, TimestepsTaken, tottime ) 
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: ParticleStepsTaken, TimestepsTaken
    REAL(KIND=dp) :: tottime

    INTEGER :: TotParticleStepsTaken, TotNoParticles

    CALL ParticleStatusCount( Particles )

    IF( ParEnv % PEs > 1 ) THEN
      TotNoParticles =  ParallelReduction( Particles % NumberOfParticles ) 
      TotParticleStepsTaken = ParallelReduction( ParticleStepsTaken) 
    ELSE
      TotNoParticles = Particles % NumberOfParticles 
      TotParticleStepsTaken =  ParticleStepsTaken
    END IF
    
    WRITE (Message,'(A,T22,I12)') 'Active particles:',TotNoParticles
    CALL Info('ParticleInformation',Message,Level=6)
    WRITE (Message,'(A,T22,ES12.2)') 'Elapsed time:',tottime
    CALL Info('ParticleInformation',Message,Level=6)
    WRITE (Message,'(A,T22,I12)') 'Time steps taken:',TimeStepsTaken
    CALL Info('ParticleInformation',Message,Level=8)
    WRITE (Message,'(A,T22,I12)') 'Particle steps taken:',TotParticleStepsTaken
    CALL Info('ParticleInformation',Message,Level=8)

  END SUBROUTINE ParticleInformation



  !---------------------------------------------------------
  !> Computes the characterestic speed for time integration.
  !> The speed may be either computed for the whole set or
  !> alternatively to just one particle.
  !---------------------------------------------------------
  FUNCTION CharacteristicSpeed( Particles, No ) RESULT ( CharSpeed )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: No
    REAL(KIND=dp) :: CharSpeed
    
    REAL(KIND=dp) :: Velo(3),Speed,SumSpeed,MaxSpeed
    INTEGER :: i,j,Cnt,NoParticles,dim,ParallelParticles
    REAL(KIND=dp), POINTER :: Velocity(:,:)
    INTEGER, POINTER :: Status(:)
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found, UseMaxSpeed, Visited = .FALSE.
    
    SAVE Visited, MaxSpeed
    
    IF(.NOT. Visited ) THEN
      Params => ListGetSolverParams()
      UseMaxSpeed = GetLogical( Params,'Characteristic Max Speed',Found)
      Visited = .TRUE.
    END IF
      
    dim = Particles % dim
    Velocity => Particles % Velocity
    
    IF( PRESENT(No)) THEN
      Velo(1:dim) = Velocity(No,1:dim)
      CharSpeed = SQRT( SUM( Velo(1:dim) ** 2 ) )
      RETURN
    END IF
    
    NoParticles =  Particles % NumberOfParticles
    Status => Particles % Status
    CharSpeed = 0.0_dp
    Velo = 0.0_dp
    Cnt = 0
    
    ! Compute characteristic speed for square since it avoids taking the sqrt
    DO i=1,NoParticles
      IF( Status(i) >= PARTICLE_LOST ) CYCLE
      IF( Status(i) < PARTICLE_INITIATED ) CYCLE
      
      Cnt = Cnt + 1
      Velo(1:dim) = Velocity(i,1:dim)
	
      Speed = SUM( Velo(1:dim) ** 2 ) 
      IF( UseMaxSpeed ) THEN
        MaxSpeed = MAX( MaxSpeed, Speed ) 
      ELSE
        SumSpeed = SumSpeed + Speed
      END IF
    END DO
    
    IF( Cnt == 0 ) RETURN
    
    IF( UseMaxSpeed ) THEN
      CharSpeed = ParallelReduction( MaxSpeed, 2 )
    ELSE
      ParallelParticles = ParallelReduction( Cnt ) 
      CharSpeed = ParallelReduction( SumSpeed ) / ParallelParticles
    END IF
    CharSpeed = SQRT( CharSpeed ) 
    CharSpeed = MAX( Charspeed, TINY( CharSpeed ) )
    
    WRITE( Message,'(A,E13.6)') 'Speed for timestep control:',CharSpeed
    CALL Info('CharacteristicSpeed',Message,Level=12)
    
  END FUNCTION CharacteristicSpeed



  !---------------------------------------------------------
  !> Computes the characterestic time spent in an element
  !> Currently computed just for one element as computing the 
  !> size of element is a timeconsuming operation.
  !---------------------------------------------------------
  FUNCTION CharacteristicElementSize( Particles, No ) RESULT ( ElementSize )
    
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: CharTime
    INTEGER, OPTIONAL :: No
    
    REAL(KIND=dp) :: ElementSize, u, v, w, DetJ, h0, h
    REAL(KIND=dp) :: ElementSizeMin, ElementSizeMax, ElementSizeAve
    REAL(KIND=dp), POINTER :: Basis(:), SizeValues(:)
    LOGICAL :: Stat, ConstantDt, UseMinSize, Found, Visited = .FALSE.
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: i, t, n, NoElems, dim
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Variable_t), POINTER :: ElementSizeVar
    TYPE(ValueList_t), POINTER :: Params
    
    SAVE Visited, Mesh, dim, Nodes, Basis, h0, SizeValues

    ConstantDt = .NOT. PRESENT( No ) 
    
    
    IF( Visited ) THEN
      IF( ConstantDt ) THEN
        ElementSize = h0       
      ELSE
        i = Particles % ElementIndex(No)        
        IF( i > 0 ) THEN
          ElementSize = SizeValues(i)
        ELSE
          ElementSize = 0._dp
        END IF
      END IF
      RETURN
    END IF
    
    Mesh => GetMesh()
    n = Mesh % MaxElementNodes
    dim = Mesh % MeshDim
    ALLOCATE( Basis(n) )
    NoElems = Mesh % NumberOfBulkElements    
 
    IF( .NOT. ConstantDt ) THEN
      ! Use existing variable only if it is of correct size!
      ElementSizeVar => VariableGet( Mesh % Variables,'Element Size' )
      NULLIFY( SizeValues )

      IF( ASSOCIATED( ElementSizeVar ) ) THEN
        SizeValues => ElementSizeVar % Values
        IF( ASSOCIATED( SizeValues ) ) THEN
          IF( SIZE( SizeValues ) /= NoElems ) NULLIFY( SizeValues )
        END IF
      END IF
      IF( .NOT. ASSOCIATED( SizeValues ) ) THEN
        ALLOCATE( SizeValues( NoElems ) )
      END IF
      SizeValues = 0.0_dp
    END IF
      
    ElementSizeMin = HUGE( ElementSizeMin ) 
    ElementSizeMax = 0.0_dp
    ElementSizeAve = 0.0_dp
    NoElems = Mesh % NumberOfBulkElements    
    
    DO t=1,NoElems
      Element => Mesh % Elements(t)      
      CALL GetElementNodes( Nodes, Element ) 
      n = GetElementNOFNodes()

      IP = GaussPoints( Element )
      u = SUM( IP % u ) / IP % n 
      v = SUM( IP % v ) / IP % n 
      w = SUM( IP % w ) / IP % n 

      stat = ElementInfo( Element, Nodes, U, V, W, detJ, Basis )
      h = detJ ** (1.0_dp / dim ) 
      ElementSizeMin  = MIN( ElementSizeMin, h )
      ElementSizeMax  = MAX( ElementSizeMax, h )
      ElementSizeAve  = ElementSizeAve + h 

      IF( .NOT. ConstantDt ) SizeValues(t) = h     
    END DO

    ElementSizeMin = ParallelReduction( ElementSizeMin, 1 ) 
    ElementSizeMax = ParallelReduction( ElementSizeMax, 2 ) 
    ElementSizeAve = ParallelReduction( ElementSizeAve ) 
    NoElems = ParallelReduction( NoElems ) 
    
    ElementSizeAve = ElementSizeAve / NoElems 

    WRITE(Message,'(A,ES12.3)') 'Minimum element size:',ElementSizeMin
    CALL Info('CharacteristicElementSize', Message,Level=12)

    WRITE(Message,'(A,ES12.3)') 'Maximum element size:',ElementSizeMax
    CALL Info('CharacteristicElementSize', Message,Level=12)

    WRITE(Message,'(A,ES12.3)') 'Average element size:',ElementSizeAve
    CALL Info('CharacteristicElementSize', Message,Level=12)

    Params => ListGetSolverParams()
    UseMinSize = GetLogical( Params,'Characteristic Minimum Size',Found)
    IF( UseMinSize ) THEN
      h0 = ElementSizeMin
    ELSE
      h0 = ElementSizeAve
    END IF

    IF(ConstantDt) THEN
       ElementSize =  h0
    ELSE
      ElementSize =  SizeValues(No)
    END IF
          
    Visited = .TRUE.
    
  END FUNCTION CharacteristicElementSize



  !-----------------------------------------------------------------------
  !> Computes the characterestic time spent for each direction separately.
  !> Currently can only be computed for one particle at a time. 
  !-----------------------------------------------------------------------
  FUNCTION CharacteristicUnisoTime( Particles, No ) RESULT ( CharTime )
    
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: CharTime
    INTEGER, OPTIONAL :: No
    
    REAL(KIND=dp) :: Center(3),CartSize(3),Velo(3)
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: i, t, n, dim
    LOGICAL :: Visited = .FALSE.
    
    SAVE Visited, Mesh, dim, Nodes

    IF(.NOT. Visited ) THEN
      Mesh => GetMesh()
      dim = Mesh % MeshDim
      Visited = .TRUE.
    END IF

    ! Absolute of velocity in each direction
    Velo(1:dim) = ABS( Particles % Velocity(No,1:dim) )

    t = Particles % ElementIndex(No)               
    Element => Mesh % Elements(t)      
    CALL GetElementNodes( Nodes, Element ) 
    n = Element % TYPE % NumberOfNodes

    ! Center point of element
    Center(1) = SUM( Nodes % x(1:n) ) / n
    Center(2) = SUM( Nodes % y(1:n) ) / n
    Center(3) = SUM( Nodes % z(1:n) ) / n

    ! Average distance from center multiplied by two in eadh direction
    CartSize(1) = 2 * SUM( ABS( Nodes % x(1:n) - Center(1) ) ) / n
    CartSize(2) = 2 * SUM( ABS( Nodes % y(1:n) - Center(2) ) ) / n
    CartSize(3) = 2 * SUM( ABS( Nodes % z(1:n) - Center(3) ) ) / n

    CharTime = HUGE( CharTime )
    DO i=1,dim
      IF( CharTime * Velo(i) > CartSize(i) ) THEN
        CharTime = CartSize(i) / Velo(i)
      END IF      
    END DO
    
  END FUNCTION CharacteristicUnisoTime


  
  !---------------------------------------------------------
  !> Computes the characterestic time spent in an element
  !> Currently computed just for one element as computing the 
  !> size of element is a timeconsuming operation.
  !---------------------------------------------------------
  FUNCTION CharacteristicElementTime( Particles, No ) RESULT ( CharTime )
    
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: CharTime
    INTEGER, OPTIONAL :: No

    REAL(KIND=dp) :: CharSpeed, CharSize

    CharSpeed = CharacteristicSpeed( Particles, No ) 
    CharSize = CharacteristicElementSize( Particles, No ) 
    CharTime = CharSize / CharSpeed

    IF( .NOT. PRESENT( No ) ) THEN
      WRITE(Message,'(A,ES12.3)') 'Characteristic time of particle:',CharTime
      CALL Info('CharacteristicElementTime', Message,Level=10)
    END IF
      
  END FUNCTION CharacteristicElementTime
  
  
  !------------------------------------------------------------------------
  !> Finds a random point that is guaranteed within the element
  !-------------------------------------------------------------------------
  FUNCTION RandomPointInElement( Element, Nodes ) RESULT ( Coord ) 

    TYPE(Element_t) :: Element
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: Coord(3)

    REAL(KIND=dp) :: u,v,w,DetJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    INTEGER :: family,n
    LOGICAL :: Stat

    family = Element % TYPE % ElementCode / 100
    n = Element % TYPE % NumberOfNodes

    ALLOCATE( Basis(n) )


100 SELECT CASE ( family )
       
    CASE ( 2 )
      u = 2*EvenRandom() - 1.0
      
    CASE ( 3 )
      u = EvenRandom()
      v = EvenRandom()
      IF( u + v > 1.0_dp ) GOTO 100
      
    CASE ( 4 )
      u = 2*EvenRandom() - 1.0
      v = 2*EvenRandom() - 1.0
  
    CASE ( 5 )
      u = EvenRandom()
      v = EvenRandom()
      w = EvenRandom()
      IF( u + v + w > 1.0_dp ) GOTO 100

    CASE ( 8 ) 
      u = 2*EvenRandom() - 1.0
      v = 2*EvenRandom() - 1.0
      w = 2*EvenRandom() - 1.0
      
    CASE DEFAULT
      CALL Fatal('RandomPointInElement','Not implemented for elementtype')
      
    END SELECT

    Stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

    Coord(1) = SUM( Basis(1:n) * Nodes % x(1:n) )
    Coord(2) = SUM( Basis(1:n) * Nodes % y(1:n) )
    Coord(3) = SUM( Basis(1:n) * Nodes % z(1:n) )

  END FUNCTION RandomPointInElement


  !------------------------------------------------------------------------
  !> Initialize particle positions and velocities with a number of different
  !> methods, both random and uniform.
  !-------------------------------------------------------------------------
  SUBROUTINE InitializeParticles( Particles, InitParticles, AppendParticles, &
      Group, IsAdvector ) 
   
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: InitParticles
    LOGICAL, OPTIONAL :: AppendParticles
    INTEGER, OPTIONAL :: Group
    LOGICAL, OPTIONAL :: IsAdvector
    
    TYPE(ValueList_t), POINTER :: Params, BodyForce 
    TYPE(Variable_t), POINTER :: VeloVar, MaskVar, AdvVar
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t) :: Nodes
    INTEGER :: Offset, NewParticles,LastParticle,NoElements
    INTEGER :: dim, ElementIndex, body_id, bf_id
    REAL(KIND=dp), POINTER :: rWork(:,:),Coordinate(:,:), Velocity(:,:)
    REAL(KIND=dp) :: Velo(3), Coord(3), Center(3), CenterVelo(3), time0, dist
    CHARACTER(:), ALLOCATABLE :: InitMethod
    INTEGER :: i,j,k,l,n,vdofs,nonodes, InitStatus, TotParticles, No
    INTEGER, POINTER :: MaskPerm(:), InvPerm(:), NodeIndexes(:)
    LOGICAL :: Found, GotIt, GotMask, RequirePositivity, GotWeight
    REAL(KIND=dp), POINTER :: InitialValues(:,:)
    REAL(KIND=dp) :: mass,boltz,temp,coeff,eps,frac,meanval 
    REAL(KIND=dp) :: MinCoord(3), MaxCoord(3), Diam, DetJ, MinDetJ, MaxDetJ, &
        MinWeight, MaxWeight, MeanWeight, Phi
    REAL(KIND=dp), POINTER :: MaskVal(:)
    REAL(KIND=dp), ALLOCATABLE :: Weight(:)
    INTEGER :: nx,ny,nz,nmax,ix,iy,iz,ind
    LOGICAL :: CheckForSize, Parallel, SaveParticleOrigin, AdvectorMode
    LOGICAL, POINTER :: DoneParticle(:)
    CHARACTER(:), ALLOCATABLE :: VariableName
    CHARACTER(LEN=MAX_NAME_LEN) :: AdvName

    CHARACTER(*), PARAMETER :: Caller = 'InitializeParticles'

    
    SAVE Nodes


    Mesh => GetMesh()
    Params => ListGetSolverParams()
    dim = Particles % Dim
    Parallel = ( ParEnv % PEs > 1 )

    n = MAX( Mesh % MaxElementNodes, 27 )
    IF(.NOT. ASSOCIATED( Nodes % x ) ) THEN
      ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    END IF      
    
    !------------------------------------------------------------------------
    ! Initialize the timestepping strategy stuff
    !-------------------------------------------------------------------------    
    AdvectorMode = .FALSE.
    IF(PRESENT(IsAdvector)) AdvectorMode = IsAdvector
    SaveParticleOrigin = AdvectorMode
    
    InitMethod = ListGetString( Params,'Coordinate Initialization Method',gotIt ) 
    VariableName = ListGetString( Params,'Velocity Variable Name',GotIt )
    IF(.NOT. GotIt) VariableName = 'Flow Solution'
    VeloVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
    vdofs = 0
    IF(ASSOCIATED(VeloVar)) THEN
      vdofs = VeloVar % dofs
      IF(.NOT. ASSOCIATED(VeloVar % Perm)) THEN
        CALL Fatal(Caller,'Velocity variable should be a field with Permutation!')
      END IF
      IF(vdofs < dim) THEN
        CALL Fatal(Caller,'Dimension of velocity variable is too small: '//I2S(vdofs))
      END IF
    END IF
    
    VariableName = ListGetString( Params,'Advector Variable Name',GotIt)
    IF(.NOT. GotIt) VariableName = 'AdvectorData'
    AdvVar => VariableGet( Mesh % Variables, VariableName )
    IF(AdvectorMode) THEN
      GotIt = .FALSE.
      IF( ASSOCIATED( AdvVar ) ) THEN
        IF( AdvVar % TYPE == Variable_on_elements ) THEN
          InitMethod = 'advector elemental'
          CALL Info(Caller,&
              'Initializing particles on center of bulk elements',Level=10)
          GotIt = .TRUE.
        ELSE IF( AdvVar % TYPE == Variable_on_nodes_on_elements ) THEN 
          InitMethod = 'advector dg'
          CALL Info(Caller,&
              'Initializing particles on scaled dg points',Level=10)
          GotIt = .TRUE.
        ELSE IF( AdvVar % TYPE == Variable_on_gauss_points ) THEN
          InitMethod = 'advector ip'
          CALL Info(Caller,&
              'Initializing particles on gauss points',Level=10)
          GotIt = .TRUE.
        END IF
      END IF
      IF(.NOT. GotIt) THEN
        InitMethod = 'advector nodal'          
        CALL Info(Caller,'Initializing particles on each node',Level=10)
      END IF
      IF( vdofs == 0) THEN
        CALL Fatal(Caller,'Velocity variable needed to initialize advector')
      END IF
    END IF
    
    Particles % RK2 = ListGetLogical( Params,'Runge Kutta', GotIt )
    IF( Particles % RK2 .AND. ParEnv % PEs > 1 ) THEN
      CALL Warn(Caller,'> Runge Kutta < integration might not work in parallel')
    END IF
    
    Particles % DtConstant = ListGetLogical( Params,'Particle Dt Constant',GotIt )
    IF(.NOT. GotIt) Particles % DtConstant = .TRUE.
    
    !------------------------------------------------------------------------
    ! The user may use a mask to initialize the particles only at a part of the 
    ! domain, or to utilize the ordeing of the permutation vector.
    ! Create the mask before deciding on the number which may be relative
    !-------------------------------------------------------------------------  
    GotMask = .FALSE.
    VariableName = ListGetString( Params,'Initialization Condition Variable',GotIt )
    IF(GotIt) THEN
      RequirePositivity = .TRUE.
    ELSE
      VariableName = ListGetString( Params,'Initialization Mask Variable',GotIt )
      RequirePositivity = .FALSE.
    END IF

    IF(GotIt) THEN
      MaskVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
      IF( .NOT. ASSOCIATED( MaskVar ) ) THEN
        CALL Fatal(Caller,'Mask / Condition variable does not exist!')
      END IF

      MaskPerm => MaskVar % Perm
      MaskVal => MaskVar % Values

      IF(.NOT. ( ASSOCIATED( MaskPerm ) .AND. ASSOCIATED(MaskVal)) ) THEN
        CALL Warn(Caller,'Initialization variable does not exist?')
      ELSE IF( MAXVAL( MaskPerm ) == 0 ) THEN
        CALL Warn(Caller,'Initialization variable of size zero?')
        nonodes = 0
        noelements = 0
        InvPerm => NULL()
      ELSE
        GotMask = .TRUE.
        IF( InitMethod == 'nodal') THEN
          ALLOCATE( InvPerm(SIZE(MaskPerm)) )
          InvPerm = 0
          j = 0
          DO i=1,SIZE(MaskPerm)
            k = MaskPerm(i)
            IF( k == 0 ) CYCLE
            IF( RequirePositivity ) THEN
              IF( MaskVal( k ) < 0.0_dp ) CYCLE
            END IF
            j = j + 1
            InvPerm(j) = i
          END DO
          nonodes = j

          CALL Info(Caller,'Total nodes '//I2S(Mesh % NumberOfNodes)//&
              ' and masked nodes '//I2S(nonodes),Level=10)
        ELSE IF( InitMethod == 'elemental') THEN
          ALLOCATE( InvPerm( MAX( Mesh % NumberOfBulkElements, Mesh % NumberOfBoundaryElements ) ) ) 
          InvPerm = 0

          j = 0
          DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
            Element => Mesh % Elements(i)
            NodeIndexes =>  Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes

            IF( i == Mesh % NumberOfBulkElements ) THEN
              IF( j > 0 ) EXIT
            END IF

            IF( ANY( MaskPerm( NodeIndexes ) == 0 ) ) CYCLE

            IF( RequirePositivity ) THEN
              meanval = SUM( MaskVal( MaskPerm( NodeIndexes ) ) ) 
              IF( meanval < 0.0_dp ) CYCLE
            END IF

            ! If some of bulk elements have been found active
            j = j + 1
            InvPerm(j) = i

          END DO
          noelements = j
          CALL Info(Caller,'Total elements '//I2S(Mesh % NumberOfBulkElements)//&
              ' and masked elements '//I2S(noelements),Level=10)
        END IF
      END IF
    ELSE
      nonodes = Mesh % NumberOfNodes
      noelements = Mesh % NumberOfBulkElements
    END IF
    
    GotWeight = ListCheckPresentAnyBodyForce(CurrentModel, &
        'Particle Initialization Weight')
    IF( GotWeight ) THEN
      CALL Info(Caller,'Using weight when creating particles',Level=8)
      ALLOCATE( Weight( Mesh % MaxElementNodes ) )
      Weight = 0.0_dp
    END IF


    !------------------------------------------------------------------------
    ! Use a simple bounding box for initialization
    ! By default a local bounding box is used...
    !-------------------------------------------------------------------------  
    IF(LEN(InitMethod)>=3) THEN
      IF( InitMethod(1:3) == 'box') THEN
        Eps = GetCReal( Params,'Wall Particle Radius',GotIt)
        IF(.NOT. GotIt) eps = 1.0d-8
      
        MinCoord(1) = GetCReal( Params,'Min Initial Coordinate 1',GotIt) 
        IF(.NOT. GotIt) MinCoord(1) = Particles % LocalMinCoord(1) + eps
      
        MaxCoord(1) = GetCReal( Params,'Max Initial Coordinate 1',GotIt) 
        IF(.NOT. GotIt) MaxCoord(1) = Particles % LocalMaxCoord(1) - eps
      
        MinCoord(2) = GetCReal( Params,'Min Initial Coordinate 2',GotIt) 
        IF(.NOT. GotIt) MinCoord(2) = Particles % LocalMinCoord(2) + eps
      
        MaxCoord(2) = GetCReal( Params,'Max Initial Coordinate 2',GotIt) 
        IF(.NOT. GotIt) MaxCoord(2) = Particles % LocalMaxCoord(2) - eps
      
        MinCoord(3) = GetCReal( Params,'Min Initial Coordinate 3',GotIt) 
        IF(.NOT. GotIt) MinCoord(3) = Particles % LocalMinCoord(3) 
      
        MaxCoord(3) = GetCReal( Params,'Max Initial Coordinate 3',GotIt) 
        IF(.NOT. GotIt) MaxCoord(3) = Particles % LocalMaxCoord(3) - eps
      END IF
    END IF
    
    IF( InitMethod == 'box random cubic' .OR. InitMethod == 'box uniform cubic') THEN
      Diam = 2 * GetCReal( Params,'Particle Cell Radius',GotIt)
      IF(.NOT. GotIt ) THEN
        Diam = GetCReal( Params,'Particle Cell Size',GotIt)
      END IF
      IF(.NOT. GotIt ) THEN
        Diam = 2 * GetCReal( Params,'Particle Radius',GotIt)
      END IF
      IF(.NOT. GotIt ) THEN
        CALL Fatal(Caller,'Size of unit cell not given')
      END IF
      
      nx = NINT ( ( MaxCoord(1) - MinCoord(1) ) / Diam )
      ny = NINT( ( MaxCoord(2) - MinCoord(2) ) / Diam )
      IF( dim == 3 ) THEN
        nz = NINT( ( MaxCoord(3) - MinCoord(3) ) / Diam )
      ELSE
        nz = 1
      END IF
    END IF

    !------------------------------------------------------------------------
    ! Now decide on the number of particles.
    !-------------------------------------------------------------------------  
    IF( AdvectorMode ) THEN
      IF( ASSOCIATED( AdvVar ) ) THEN
        NewParticles = SIZE( AdvVar % Values )
      ELSE
        NewParticles = Mesh % NumberOfNodes
      END IF
      CALL Info(Caller,'Using ParticleData variable of size '&
          //I2S(NewParticles)//' to define particles!')
    ELSE IF( PRESENT( InitParticles ) ) THEN
      NewParticles = InitParticles
    ELSE
      IF( InitMethod == 'box uniform cubic') THEN
        NewParticles = nx * ny * nz
      ELSE
        NewParticles = GetInteger( Params,'Number of Particles',GotIt) 
      END IF
      IF(.NOT. GotIt ) THEN
        frac = GetCReal( Params,'Particle Node Fraction',GotIt)      
        IF( GotIt ) THEN
          NewParticles = NINT( frac * nonodes )
        ELSE
          frac = GetCReal( Params,'Particle Element Fraction',GotIt)
          IF( GotIt ) THEN
            NewParticles = NINT( frac * noelements )
          ELSE
            frac = GetCReal( Params,'Particle Cell Fraction',GotIt)
            IF( GotIt ) THEN
              NewParticles = NINT( frac * nx * ny * nz )
            ELSE
              CALL Fatal(Caller,'Could not determine the number of new particles!')
            END IF
          END IF
        END IF
      END IF
    END IF
    
    IF( ParEnv% PEs == 1 ) THEN
      TotParticles = NewParticles
    ELSE
      TotParticles = ParallelReduction( NewParticles ) 
    END IF
    
    IF( TotParticles == 0 ) THEN
      CALL Fatal(Caller,'No Particles to Initialize')
    ELSE
      WRITE( Message,'(A,I0)') 'Number of Particles: ',TotParticles
      CALL Info(Caller,Message,Level=6)
    END IF

    !------------------------------------------------------------------------
    ! If there are no particles in this partition, nothing to do
    !------------------------------------------------------------------------- 
    IF( NewParticles == 0 ) RETURN
        
    !------------------------------------------------------------------------
    ! Interval of particles
    !-------------------------------------------------------------------------  
    IF( PRESENT( AppendParticles ) ) THEN
      Offset = Particles % NumberOfParticles
    ELSE
      Offset = 0
    END IF
    LastParticle = Offset + NewParticles
    
    
    !------------------------------------------------------------------------
    ! Allocate particles
    !-------------------------------------------------------------------------    
    CALL AllocateParticles( Particles, LastParticle )

    IF( SaveParticleOrigin ) THEN   
      IF(.NOT. ASSOCIATED( Particles % NodeIndex ) ) THEN
        ALLOCATE( Particles % NodeIndex( NewParticles ) )
        Particles % NodeIndex = 0
      END IF
       
      IF( Parallel ) THEN
        IF(.NOT. ASSOCIATED( Particles % Partition ) ) THEN
          ALLOCATE( Particles % Partition( NewParticles ) )
        END IF
        Particles % Partition = ParEnv % MyPe + 1
      END IF
    END IF
        
    
    IF( Particles % NumberOfGroups > 0 ) THEN
      IF( .NOT. PRESENT( Group ) ) THEN
        CALL Fatal(Caller,'Group used inconsistently!')
      END IF
      Particles % Group(Offset+1:LastParticle) = Group
    END IF
      
    Particles % NumberOfParticles = LastParticle
    
    Velocity => Particles % Velocity
    Coordinate => Particles % Coordinate
    
    
    SELECT CASE ( InitMethod ) 
      
    CASE ('nodal ordered')
      CALL Info(Caller,'Initializing particles evenly among nodes',Level=10)

      Particles % NumberOfParticles = NewParticles
      DO i=1,NewParticles
        k = Offset + i
        j = (nonodes-1)*(i-1)/(NewParticles-1)+1
        IF( GotMask ) j = InvPerm(j)
        Coordinate(k,1) = Mesh % Nodes % x(j)
        Coordinate(k,2) = Mesh % Nodes % y(j)
        IF( dim == 3 ) Coordinate(k,3) = Mesh % Nodes % z(j)
      END DO
      
    CASE ('elemental random')
      CALL Info(Caller,&
          'Initializing particles randomly within elements',Level=10)

      Particles % NumberOfParticles = NewParticles

      MaxDetJ = 0.0_dp
      MinDetJ = HUGE( MinDetJ )
      MaxWeight = -HUGE( MaxWeight ) 
      MinWeight = HUGE( MinWeight )

      DO i = 1, NoElements

        j = i
        IF( GotMask ) j = InvPerm(j)

        Element => Mesh % Elements(j)
        NodeIndexes =>  Element % NodeIndexes
        n = Element % TYPE % NumberOfNodes

        ! If weight is used see that we have a weight, and that it is positive
        IF( GotWeight ) THEN
          IF( j > Mesh % NumberOfBulkElements ) CYCLE

          body_id = Element % BodyId
          bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values,&
              'Body Force',minv=1)
          BodyForce => CurrentModel % BodyForces(bf_id) % Values
          Weight(1:n) = ListGetReal( BodyForce,'Particle Initialization Weight',&
              n, NodeIndexes, GotIt)
          IF(.NOT. GotIt) CYCLE

          MeanWeight = SUM( Weight(1:n) ) / n
          MaxWeight = MAX( MaxWeight, MeanWeight )
          MinWeight = MIN( MinWeight, MeanWeight ) 
          IF( MeanWeight <= 0.0_dp ) CYCLE
        END IF

        ! Compute the size of the element if this is an active element.
        Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

        DetJ = ElementSize( Element, Nodes ) 
        MaxDetJ = MAX( MaxDetJ, DetJ ) 
        MinDetJ = MIN( MinDetJ, DetJ )
      END DO      
      
      WRITE( Message,'(A,ES12.3)') 'Maximum size of elements:',MaxDetJ
      CALL Info(Caller,Message,Level=8)
      WRITE( Message,'(A,ES12.3)') 'Minimum size of elements:',MinDetJ
      CALL Info(Caller,Message,Level=8)
      IF( GotWeight ) THEN
        WRITE( Message,'(A,ES12.3)') 'Maximum weight in elements:',MaxWeight
        CALL Info(Caller,Message,Level=8)
        WRITE( Message,'(A,ES12.3)') 'Minimum weight in elements:',MinWeight
        CALL Info(Caller,Message,Level=8)
      END IF

      IF( MaxWeight < 0.0 ) THEN
        CALL Info(Caller,'No positive weight!')
        RETURN
      END IF

      ! If all elements are of same size then no need to check for size
      CheckForSize = ( MinDetJ < (1-EPSILON(MaxDetJ))*MaxDetJ )
      IF( GotWeight ) THEN
        CheckForSize = CheckForSize .OR. &
            ( MinWeight < (1-EPSILON(MaxWeight))*MaxWeight )
      END IF

      i = 0      
      DO WHILE(.TRUE.) 

        j = CEILING( NoElements * EvenRandom() )
        IF( GotMask ) j = InvPerm(j)
                
        Element => Mesh % Elements(j)
        NodeIndexes =>  Element % NodeIndexes
        n = Element % TYPE % NumberOfNodes

        Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        IF( CheckForSize ) THEN
          DetJ = ElementSize( Element, Nodes ) 

          ! The weight could be computed really using the integration point
          ! Here we assumes constant weight within the whole element. 
          IF( GotWeight ) THEN
            IF( j > Mesh % NumberOfBulkElements ) CYCLE

            body_id = Element % BodyId
            bf_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values,&
                'Body Force',minv=1)
            BodyForce => CurrentModel % BodyForces(bf_id) % Values
            Weight(1:n) = ListGetReal( BodyForce,'Particle Initialization Weight',&
                n, NodeIndexes, GotIt)
            IF(.NOT. GotIt ) CYCLE

            MeanWeight = SUM( Weight(1:n) ) / n
            IF( MeanWeight <= 0.0_dp ) CYCLE

            ! Do importance samping for the particles
            IF( EvenRandom() * MaxDetJ * MaxWeight > DetJ * MeanWeight ) CYCLE
          ELSE          
            IF( EvenRandom() * MaxDetJ > DetJ ) CYCLE
          END IF
        END IF

        ! Create a random particle within the element
        Coord = RandomPointInElement( Element, Nodes )

        i = i + 1
        k = Offset + i
        Coordinate(k,1:dim) = Coord(1:dim)

        ! Only a bulk element may own a particle
        IF( j <= Mesh % NumberOfBulkElements ) THEN
          Particles % ElementIndex(k) = j
        END IF

        IF( i == NewParticles ) EXIT
      END DO

    CASE ('elemental ordered')
      CALL Info(Caller,&
          'Initializing particles evenly among elements',Level=10)

      NewParticles = MIN(NoElements,NewParticles)
      Particles % NumberOfParticles = NewParticles

      DO i=1,NewParticles
        k = Offset + i
        j = (NoElements-1)*(i-1)/(NewParticles-1)+1
        IF( GotMask ) j = InvPerm(j)
                
        IF( j > Mesh % NumberOfBulkElements ) THEN
          PRINT *,'j too large',j,i,k,(NoElements-1)*(i-1)/(NewParticles-1)+1
        END IF
        
        Element => Mesh % Elements(j)
        NodeIndexes =>  Element % NodeIndexes
        n = Element % TYPE % NumberOfNodes
        Coord(1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n
        Coord(2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n
        IF( dim == 3 ) Coord(3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n

        Coordinate(k,1:dim) = Coord(1:dim)
        
        ! Only a bulk element may own a particle
        IF( j <= Mesh % NumberOfBulkElements ) THEN
          Particles % ElementIndex(i) = j
        END IF
      END DO
      
    CASE ('advector nodal')
      DO i=1,Mesh % NumberOfNodes
        j = i
        IF( ASSOCIATED( AdvVar ) ) THEN
          j = AdvVar % Perm(i)
	  IF(j==0) CYCLE
        END IF
        Coordinate(j,1) = Mesh % Nodes % x(i)
        Coordinate(j,2) = Mesh % Nodes % y(i)
        IF( dim == 3 ) Coordinate(j,3) = Mesh % Nodes % z(i)

        IF(vdofs>0) THEN
          DO k=1,dim	
            Velocity(j,k) = VeloVar % Values(vdofs*(VeloVar % Perm(i)-1)+k)
          END DO
        END IF
      END DO

      IF( SaveParticleOrigin ) THEN
        DO i=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(i)
          NodeIndexes =>  Element % NodeIndexes
          n = Element % TYPE % NumberOfNodes
          DO j=1,n
            IF( ASSOCIATED( AdvVar ) ) THEN
              k = AdvVar % Perm(NodeIndexes(j))
            ELSE
              k = NodeIndexes(j)
            END IF
            IF(k==0) CYCLE
            Particles % ElementIndex(k) = i
            Particles % NodeIndex(k) = NodeIndexes(j)
          END DO
        END DO
      END IF
      
    CASE ('advector elemental')
      DO i=1,NoElements                       
        No = AdvVar % Perm( i )
        IF( No == 0 ) CYCLE
        
        Element => Mesh % Elements(i)
        NodeIndexes =>  Element % NodeIndexes
        n = Element % TYPE % NumberOfNodes
        
        Center(1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n
        Center(2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n
        IF( dim == 3 ) Center(3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n
        Coordinate(No,1:dim) = Center(1:dim)        

        IF(vdofs>0) THEN
          DO j=1,dim
            CenterVelo(j) = SUM( VeloVar % Values(vdofs*(VeloVar % Perm(NodeIndexes)-1)+j) ) / n
          END DO
          Velocity(No,1:dim) = CenterVelo(1:dim) 
        END IF
                  
        Particles % ElementIndex(No) = i
        IF( SaveParticleOrigin ) THEN
          Particles % NodeIndex(No) = No
        END IF
      END DO

    CASE ('advector ip','ip')

      BLOCK
        TYPE(GaussIntegrationPoints_t) :: IP
        REAL(KIND=dp) :: detJ, Basis(27)
        INTEGER :: m
        LOGICAL :: stat
        LOGICAL :: Debug
        
        Debug = .FALSE.
        
        DO i=1,NoElements                       
          No = AdvVar % Perm( i )
          IF( No == 0 ) CYCLE
          
          Element => Mesh % Elements(i)
          NodeIndexes =>  Element % NodeIndexes
          n = Element % TYPE % NumberOfNodes
          
          Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          Nodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
                    
          IF( debug ) THEN
            PRINT *,'x:',nodes % x(1:n) 
            PRINT *,'y:',nodes % y(1:n) 
            PRINT *,'z:',nodes % z(1:n) 
          END IF
          
          m = AdvVar % Perm(i+1) - AdvVar % Perm(i) 
          IF( m == 0 ) CYCLE
          
          IP = GaussPoints(Element, m )
          
          DO j = 1, IP % n 
            stat = ElementInfo( Element, Nodes, IP % v(j), IP % u(j), IP % w(j), detJ, Basis ) 
            No = AdvVar % Perm(i) + j
            
            Coord(1) = SUM( Basis(1:n) * Nodes % x(1:n) )
            Coord(2) = SUM( Basis(1:n) * Nodes % y(1:n) )
            IF( dim == 3 ) Coord(3) = SUM( Basis(1:n) * Nodes % z(1:n) )
            Coordinate(No,1:dim) = Coord(1:dim)

            IF( vdofs > 0 ) THEN
              DO l=1,dim
                Velo(l) = SUM( Basis(1:n) * &
                    VeloVar % Values(vdofs*(VeloVar % Perm(NodeIndexes)-1)+l ) )
              END DO
              Velocity(No,1:dim) = Velo(1:dim)
            END IF
              
            IF( Debug ) THEN
              PRINT *,'j:',j,m,No,Coord(1:dim),Velo(1:dim)
            END IF
                            
            Particles % ElementIndex(No) = i
            IF( SaveParticleOrigin ) THEN
              Particles % NodeIndex(No) = No                       
            END IF
          END DO
        END DO
      END BLOCK
      
    CASE ( 'advector dg','dg')      

      BLOCK
        REAL(KIND=dp) :: DgScale
        LOGICAL :: GotScale

        DGScale = ListGetCReal( Params,'DG Nodes Scale',GotScale )
        IF(.NOT. GotScale ) DgScale = 1.0 / SQRT( 3.0_dp ) 
        GotScale = ( ABS( DGScale - 1.0_dp ) > TINY( DgScale ) )

        DO i=1,NoElements                       
          Element => Mesh % Elements(i)
          NodeIndexes =>  Element % NodeIndexes
          n = Element % TYPE % NumberOfNodes

          Element => Mesh % Elements(i)
          NodeIndexes =>  Element % NodeIndexes
          n = Element % TYPE % NumberOfNodes

          IF( GotScale ) THEN
            Center(1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n
            Center(2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n
            IF( dim == 3 ) Center(3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n
            
            IF( vdofs > 0 ) THEN
              DO j=1,dim
                CenterVelo(j) = SUM( VeloVar % Values(vdofs*(VeloVar % Perm(NodeIndexes)-1)+j) ) / n
              END DO
            END IF
          END IF
            
          DO j = 1, n
            No = AdvVar % Perm( Element % DgIndexes(j) )
            IF( No == 0 ) CYCLE
            k = NodeIndexes(j)
            Coord(1) = Mesh % Nodes % x(k)
            Coord(2) = Mesh % Nodes % y(k)
            IF( dim == 3 ) Coord(3) = Mesh % Nodes % z(k)            
            
            IF( GotScale ) THEN
              Coord(1:dim) = Center(1:dim) + ( Coord(1:dim) - Center(1:dim) ) * DgScale 
            END IF            
            Coordinate(No,1:dim) = Coord(1:dim) 
            
            IF( vdofs > 0 ) THEN
              DO l=1,dim
                Velo(l) = VeloVar % Values(vdofs*(VeloVar % Perm(k)-1)+l) 
              END DO              
              IF( GotScale ) THEN
                Velo(1:dim) = CenterVelo(1:dim) + (Velo(1:dim) - CenterVelo(1:dim)) * DgScale
              END IF              
              Velocity(No,1:dim) = Velo(1:dim) 
            END IF
              
            Particles % ElementIndex(No) = i
            
            IF( SaveParticleOrigin ) THEN
              Particles % NodeIndex(No) = No                       
            END IF
          END DO
        END DO
      END BLOCK
                  
    CASE ('sphere random')
      CALL Info(Caller,&
          'Initializing particles randomly within a sphere',Level=10)

      Diam = GetCReal( Params,'Initial Sphere Radius')
      rWork => ListGetConstRealArray( Params,'Initial Sphere Center')
      IF ( ASSOCIATED(rwork) ) THEN
        Center = rWork(1:3,1)
      ELSE
        Center = 0.0_dp
      END IF
      
      i = 0
      DO WHILE (.TRUE.) 
        DO j=1,dim
          Coord(j) = Diam*(2*EvenRandom()-1)
        END DO
        ! Is the point within sphere (or circle in 2d)
        IF( SUM( Coord(1:dim)**2 ) > Diam*Diam ) CYCLE
        
        i = i + 1
        k = Offset + i      
        Coordinate(k,:) = Center + Coord(1:dim)
        IF( i == NewParticles ) EXIT
      END DO
      
    CASE ('box random')
      DO i=1,NewParticles
        k = Offset + i      
        DO j=1,dim
          Coord(j) = MinCoord(j) + (MaxCoord(j)-MinCoord(j)) * EvenRandom()
        END DO
        Coordinate(k,:) = Coord(1:dim)
      END DO
      
    CASE ('box random cubic')
      CALL Info(Caller,&
          'Initializing particles randomly in a grid',Level=10)

      nmax = nx * ny * nz
      IF( nmax < NewParticles ) THEN
        CALL Fatal(Caller,'More particles than places in unit cell')
      END IF
      
      ALLOCATE( DoneParticle(nx*ny*nz) )
      
      IF( NewParticles == nmax ) THEN
        ! if the list is full just set all true
        DoneParticle = .TRUE.
      ELSE IF( NewParticles < nmax / 2 ) THEN
        ! If there are few particles start from an empty list and count upwards
        DoneParticle = .FALSE.
        i =  0
        DO WHILE(.TRUE.) 
          ind = NINT( NewParticles * EvenRandom() + 0.5 )
          IF( .NOT. DoneParticle(i) ) THEN
            DoneParticle(ind) = .TRUE.
            i = i + 1
            IF( i == NewParticles ) EXIT
          END IF
        END DO
      ELSE    
        ! if there are many particles start from a full list and count downwards
        DoneParticle = .TRUE.
        i = nmax
        DO WHILE(.TRUE.) 
          ind = NINT( NewParticles * EvenRandom() + 0.5 )
          IF( DoneParticle(i) ) THEN
            DoneParticle(ind) = .FALSE.
            i = i - 1
            IF( i == NewParticles ) EXIT
          END IF
        END DO
      END IF
      
      ! set the coordinates 
      i = 0
      DO ix = 1, nx
        DO iy = 1, ny
          DO iz = 1, nz
            ind = nx*ny*(iz-1) + nx*(iy-1) + ix
            IF( DoneParticle(ind) ) THEN
              i = i + 1
              k = Offset + i
              Coordinate(k,1) = MinCoord(1) + ( 1.0_dp*ix - 0.5) * Diam 
              Coordinate(k,2) = MinCoord(2) + ( 1.0_dp*iy - 0.5) * Diam 
              IF( dim == 3 ) THEN
                Coordinate(k,3) = MinCoord(3) + ( 1.0_dp*iz - 0.5) * Diam 
              END IF
            END IF
          END DO
        END DO
      END DO
      DEALLOCATE( DoneParticle ) 
      
    CASE ('box uniform cubic')
      CALL Info(Caller,&
          'Initializing particles in a grid',Level=10)

      nmax = nx * ny * nz
      IF( nmax /= NewParticles ) THEN
        CALL Fatal(Caller,'Wrong number of particles')
      END IF
      
      ! set the coordinates 
      i = 0
      DO ix = 1, nx
        DO iy = 1, ny
          DO iz = 1, nz
            ind = nx*ny*(iz-1) + nx*(iy-1) + ix
            i = i + 1
            k = Offset + i
            Coordinate(k,1) = MinCoord(1) + ( 1.0_dp*ix - 0.5) * Diam 
            Coordinate(k,2) = MinCoord(2) + ( 1.0_dp*iy - 0.5) * Diam 
            IF( dim == 3 ) THEN
              Coordinate(k,3) = MinCoord(3) + ( 1.0_dp*iz - 0.5) * Diam 
            END IF
          END DO
        END DO
      END DO
     
    CASE DEFAULT       
      CALL Info(Caller,'Initializing particles using given coordinates',Level=10)

      InitialValues => ListGetConstRealArray(Params,'Initial Coordinate',gotIt)    
      IF(gotIt) THEN
        IF( SIZE(InitialValues,2) /= dim ) THEN
          CALL Fatal('ParticleTracker','Wrong dimension in > Initial Coordinate <')
        ELSE IF( SIZE(InitialValues,1) == 1 ) THEN
          DO i=1,NewParticles
            k = offset + i
            Coordinate(k,1:dim) = InitialValues(1,1:dim)
          END DO
        ELSE IF( SIZE(InitialValues,1) /= NewParticles ) THEN
          CALL Fatal('ParticleTracker','Wrong number of particles in > Initial Coordinate <')
        ELSE
          DO i=1,NewParticles
            k = Offset + i
            Coordinate(k,1:dim) = InitialValues(i,1:dim)
          END DO
        END IF
      ELSE
        CALL Fatal('ParticleTracker','> Initial Coordinate < not given')
      END IF
    END SELECT


    IF( GotMask ) THEN
      IF ( ASSOCIATED(InvPerm) ) DEALLOCATE( InvPerm ) 
    END IF

    !------------------------------------------------------------------------
    ! Velocities may be initialized using a given list, or obtaining them
    ! from random even or maxwell boltzmann distributions. These are additive to 
    ! allow bulk velocities with the random one.
    !-------------------------------------------------------------------------
    IF( .NOT. AdvectorMode ) THEN
      InitialValues => ListGetConstRealArray(Params,'Initial Velocity',gotIt)
      IF(gotIt) THEN
        IF( SIZE(InitialValues,2) /= DIM ) THEN
          CALL Fatal('ParticleTracker','Wrong dimension in Initial Velocity')
        ELSE IF( SIZE(InitialValues,1) == 1 ) THEN
          DO i=1,NewParticles
            k = Offset + i
            Velocity(k,1:dim) = InitialValues(1,1:DIM)
          END DO
        ELSE IF( SIZE(InitialValues,1) /= NewParticles ) THEN
          CALL Fatal('ParticleTracker','Wrong number of particles in Initial Velocity')
        ELSE
          DO i=1,NewParticles
            k = Offset + i
            Velocity(k,1:dim) = InitialValues(i,1:dim)
          END DO
        END IF
      END IF

      InitMethod = ListGetString( Params,'Velocity Initialization Method',gotIt ) 
      IF( GotIt ) THEN
        coeff = ListGetCReal( Params,'Initial Velocity Amplitude',GotIt)

        SELECT CASE ( InitMethod ) 

        CASE ('nodal velocity')
          CALL Info(Caller,&
              'Initializing velocities from the corresponding nodal velocity',Level=10)
          IF( vdofs == 0 ) THEN
            CALL Fatal(Caller,'Velocity variable needed for method >nodal velocity<')
          END IF

          DO i=1,NewParticles
            k = Offset + i
            l = Particles % NodeIndex(i)
            l = VeloVar % Perm(l)
            DO j=1,dim
              Velocity(k,j) = VeloVar % Values(vdofs*(l-1)+j)
            END DO
          END DO

        CASE ('elemental velocity') 
          CALL Info(Caller,&
              'Initializing velocities from the corresponding elemental velocity',Level=10)      
          IF( vdofs == 0 ) THEN
            CALL Fatal(Caller,'Velocity variable needed for method >elemental velocity<')
          END IF

          DO i=1,NewParticles
            k = Offset + i
            l = Particles % NodeIndex(i) ! now an elemental index
            NodeIndexes => Mesh % Elements(l) % NodeIndexes
            DO j=1,dim
              Velocity(k,j) = SUM( VeloVar % Values(vdofs*(VeloVar % Perm(NodeIndexes)-1)+j) ) / &
                  SIZE( NodeIndexes )
            END DO
          END DO

        CASE ('thermal random')  
          CALL Info(Caller,&
              'Initializing velocities from a thermal distribution',Level=10)

          IF(.NOT. GotIt) THEN
            mass = ListGetConstReal( Params,'Particle Mass')
            temp = ListGetConstReal( Params,'Particle Temperature')
            boltz = ListGetConstReal( CurrentModel % Constants,'Boltzmann constant')
            coeff = SQRT(boltz * temp / mass )
          END IF

          DO i=1,NewParticles
            k = Offset + i
            DO j=1,dim
              Velo(j) = NormalRandom()
            END DO
            Velocity(k,:) = Velocity(k,:) + coeff * Velo(1:dim)
          END DO

        CASE ('even random')
          CALL Info(Caller,'Initializing velocities from a even distribution',Level=10)

          DO i=1,NewParticles
            k = Offset + i
            DO j=1,dim
              Velo(j) = (2*EvenRandom()-1)
            END DO
            Velocity(k,:) = Velocity(k,:) + coeff * Velo(1:dim)
          END DO

        CASE ('constant random')
          CALL Info(Caller,'Initializing constant velocities in random direction',Level=10)

          DO i=1,NewParticles
            k = Offset + i
            DO j=1,dim
              Velo(j) =  (2*EvenRandom()-1)
            END DO
            Velo(1:dim) = Velo(1:dim) / SQRT(SUM(Velo(1:dim)**2))
            Velocity(k,:) = Velocity(k,:) + coeff * Velo(1:dim)
          END DO

        CASE ('constant 2d')
          CALL Info(Caller,'Initializing constant velocities evenly to space',Level=10)

          DO i=1,NewParticles
            k = Offset + i
            Phi = 2.0_dp * PI * i / NewParticles
            Velo(1) = coeff * COS( Phi )
            Velo(2) = coeff * SIN( Phi )
            Velocity(k,:) = Velocity(k,:) + Velo(1:dim)
          END DO
          
        CASE DEFAULT
          CALL Fatal(Caller,'Unknown velocity initialization method: '//TRIM(InitMethod))
          
        END SELECT
      END IF
        

      ! There may be a timestep related to initial velocity,
      ! which may be used to have the initial status developed
      ! from the initial coordinates.
      !-------------------------------------------------------
      time0 = ListGetCReal(Params,'Initial Velocity Time',gotIt)
      IF( GotIt ) THEN
        DO i=1,NewParticles
          k = Offset + i
          Coord(1:dim) = time0 * Velocity(k,:)
          Coordinate(k,:) = Coordinate(k,:) + Coord(1:dim)

          dist = SQRT( SUM( Coord(1:dim)**2 ) )
          !        Particles % Distance(k) = dist
        END DO
      END IF
    END IF

      
    ! Initialize coordinate with octree if requested
    !-------------------------------------------------------
    IF( ListGetLogical(Params,'Initial Coordinate Search',gotIt) ) THEN
      Coord = 0.0_dp
      DO i=1,NewParticles
        k = Offset + i
        ElementIndex = Particles % ElementIndex(k)
        IF( ElementIndex > 0 ) CYCLE       
        Coord(1:dim) = Coordinate(k,:) 
        CALL LocateParticleInMeshOctree( ElementIndex, Coord )
        Particles % ElementIndex(k) = ElementIndex
      END DO
    END IF
    
    !------------------------------------------------------
    ! The initial status of particles is different if using 
    ! gradual release strategy. 
    !-------------------------------------------------------
    IF( ListCheckPresent( Params,'Particle Release Number') .OR. &
      ListCheckPresent( Params,'Particle Release Fraction') ) THEN
      InitStatus = PARTICLE_WAITING
    ELSE
      InitStatus = PARTICLE_INITIATED
    END IF
  
    DO i=1,NewParticles
      k = Offset + i
      Particles % Status(k) = InitStatus
    END DO
    
    Particles % PrevCoordinate = Coordinate
    IF( ASSOCIATED( Particles % PrevVelocity ) ) THEN
      Particles % PrevVelocity = Velocity
    END IF

    ! Finally, add additional particle variables that are used for path integrals
    IF( ListCheckPresentAnyBodyForce( CurrentModel,&
        'Particle Distance Integral Source') ) THEN
      CALL ParticleVariableCreate( Particles,'Particle Distance Integral' )
      CALL ParticleVariableInitialize( Particles, Mesh,'Particle Distance Integral')
    END IF
    IF( ListCheckPresentAnyBodyForce( CurrentModel,&
        'Particle Time Integral Source') ) THEN
      CALL ParticleVariableCreate( Particles,'Particle Time Integral' )
      CALL ParticleVariableInitialize( Particles, Mesh,'Particle Time Integral')
    END IF
    
  END SUBROUTINE InitializeParticles

 

  !---------------------------------------------------------------------------
  !> This subroutine finds the possible intersection between elementfaces 
  !> and a line segment defined by two coordinates.
  !---------------------------------------------------------------------------
  SUBROUTINE SegmentElementIntersection(Mesh,BulkElement,&
      Rinit,Rfin,MinLambda,FaceElement)
    !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER   :: BulkElement
    REAL(KIND=dp) :: Rinit(3), Rfin(3), MinLambda
    TYPE(Element_t), POINTER :: FaceElement
    
    TYPE(Element_t), POINTER :: FaceElement2
    TYPE(Element_t), POINTER   :: BoundaryElement
    TYPE(Nodes_t), SAVE :: BoundaryNodes
    REAL(KIND=dp) :: Lambda(6), PosEps, NegEps, Dist
    INTEGER :: ElemDim, NoFaces, FaceIndex(6), i,j,n
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Success, AtBoundary, AtFace, Visited = .FALSE.
        
   
    PosEps = 1.0e-10
    NegEps = -1.0e-7

    ElemDim = BulkElement % TYPE % DIMENSION

    IF( ElemDim == 3 ) THEN
      NoFaces = BulkElement % TYPE % NumberOfFaces 
    ELSE 
      NoFaces = BulkElement % TYPE % NumberOfEdges 
    END IF

    DO i=1, NoFaces  
      IF( ElemDim == 3 ) THEN
        j = BulkElement % FaceIndexes(i)
        BoundaryElement => Mesh % Faces( j )
      ELSE
        j = BulkElement % EdgeIndexes(i)
        BoundaryElement => Mesh % Edges(j)
      END IF

      CALL GetElementNodes(BoundaryNodes,BoundaryElement)
        
      Lambda(i) = LineFaceIntersection(BoundaryElement,BoundaryNodes,&
          Rinit,Rfin) 
      FaceIndex(i) = j
    END DO

    ! Sort the lambdas so that the best intersection is easily found
    !-----------------------------------------------------------------------
    CALL SortR( NoFaces, FaceIndex, Lambda ) 

    ! Either there must be a positive lambda, or then the initial node must 
    ! already sit on the face.
    !------------------------------------------------------------------------
    j = 0
    DO i=1,NoFaces
      IF( Lambda(i) >= PosEps ) THEN
        j = i
      ELSE 
        IF( j > 0 ) THEN
	  MinLambda = Lambda( j ) 
          EXIT
        ELSE IF( Lambda(i) >= NegEps ) THEN
	  MinLambda = MAX( Lambda( i ), 0.0_dp )
          j = i
        END IF
        EXIT
      END IF
    END DO

    IF( j > 0 ) THEN
      IF( ElemDim == 3 ) THEN
        FaceElement => Mesh % Faces( FaceIndex(j) )
      ELSE
        FaceElement => Mesh % Edges( FaceIndex(j) )
      END IF
    ELSE 
      MinLambda = HUGE( MinLambda ) 
      FaceElement => NULL()

      CALL Warn('SegmentElementIntersection','Could not find any intersection')
      PRINT *,'Lambda: ',NoFaces,Lambda(1:NoFaces)
    END IF

!PRINT *,'Lambda:',Lambda(1:NoFaces)
!PRINT *,'MinLambda',MinLambda
    
  END SUBROUTINE SegmentElementIntersection
  

  !---------------------------------------------------------------------------
  !> This subroutine finds the possible intersection between elementfaces 
  !> and a line segment defined by two coordinates.
  !---------------------------------------------------------------------------
  SUBROUTINE SegmentElementIntersection2(Mesh,BulkElement,&
      Rinit,Rfin,MinLambda,FaceElement)
    !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER   :: BulkElement
    REAL(KIND=dp) :: Rinit(3), Rfin(3), MinLambda
    TYPE(Element_t), POINTER :: FaceElement
    
    TYPE(Element_t), POINTER :: FaceElement2
    TYPE(Element_t), POINTER   :: BoundaryElement
    TYPE(Nodes_t), SAVE :: BoundaryNodes
    REAL(KIND=dp) :: Lambda, PosEps, NegEps, Dist
    INTEGER :: ElemDim, NoFaces, i,j,n
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Success, AtBoundary, AtFace, Visited = .FALSE.
        
   
    PosEps = 1.0e-10
    NegEps = -1.0e-7
    
    ElemDim = BulkElement % TYPE % DIMENSION
    MinLambda = -HUGE( MinLambda ) 

    IF( ElemDim == 3 ) THEN
      NoFaces = BulkElement % TYPE % NumberOfFaces 
    ELSE
      NoFaces = BulkElement % TYPE % NumberOfEdges 
    END IF

    DO i=1, NoFaces          
     IF( ElemDim == 3 ) THEN
        j = BulkElement % FaceIndexes(i)
        BoundaryElement => Mesh % Faces( j )
      ELSE
        j = BulkElement % EdgeIndexes(i)
        BoundaryElement => Mesh % Edges(j)
      END IF
      
      CALL GetElementNodes(BoundaryNodes,BoundaryElement)
      
      Lambda = LineFaceIntersection2(BoundaryElement,BoundaryNodes,&
          Rinit,Rfin,Success)
      IF(.NOT. Success ) CYCLE 

      IF( Lambda > MinLambda ) THEN
        MinLambda = Lambda 
        FaceElement => BoundaryElement 
        IF( MinLambda > PosEps ) EXIT
      END IF
    END DO

    IF( MinLambda < NegEps ) THEN
      FaceElement => NULL()
!      CALL Warn('SegmentElementIntersection','Could not find any intersection')
!      PRINT *,'Lambda: ',NoFaces,MinLambda
    ELSE 
      MinLambda = MAX( MinLambda, 0.0_dp )
    END IF
    
  END SUBROUTINE SegmentElementIntersection2
  

  !---------------------------------------------------------------------------
  !> This subroutine tests whether a particle is within element using the
  !> consistent strategy with the above algorithm.
  !---------------------------------------------------------------------------
  FUNCTION SegmentElementInside(Mesh,BulkElement,Rfin,Debug) RESULT ( Inside )
    !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER   :: BulkElement
    REAL(KIND=dp) :: Rfin(3)
    LOGICAL :: Debug
    LOGICAL :: Inside
    !---------------------------------------------------------------------------
    REAL(KIND=dp) :: Rinit(3), MinLambda
    TYPE(Element_t), POINTER   :: BoundaryElement
    TYPE(Nodes_t), SAVE :: Nodes
    REAL(KIND=dp) :: Lambda, Eps, Lambdas(6)
    INTEGER :: ElemDim, NoFaces, i,j,n,Hits
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Success, AtBoundary, AtFace, Visited = .FALSE., Outside
    REAL(KIND=dp) :: minx,maxx,dx,mindist,dist


    Inside = .FALSE.

    MinLambda = HUGE(MinLambda)   
    Eps = EPSILON(Eps) !Eps =1.0e-12

    n = GetElementNOFNOdes(BulkElement)
    CALL GetElementNodes(Nodes,BulkElement)
    ElemDim = BulkElement % TYPE % DIMENSION  

    ! First make a quick and dirty test to save time
    !--------------------------------------------------------------------------
    DO i=1,ElemDim
      IF( i == 1 ) THEN
        maxx = MAXVAL( Nodes % x(1:n) )
        minx = MINVAL( Nodes % x(1:n) )
      ELSE IF( i == 2 ) THEN
        maxx = MAXVAL( Nodes % y(1:n) )
        minx = MINVAL( Nodes % y(1:n) )
      ELSE
        maxx = MAXVAL( Nodes % z(1:n) )
        minx = MINVAL( Nodes % z(1:n) )
      END IF

      dx = Eps * (maxx-minx) 
      Outside = ( Rfin(i) < minx - dx .OR. Rfin(i) > maxx + dx ) 

      IF( Debug ) THEN
        PRINT *,'Rough test: ',Outside,i,dx,minx,maxx,Rfin(i),MAX(minx-Rfin(i),Rfin(i)-maxx)
      END IF

      IF( Outside ) RETURN
    END DO


    ! Then the more laborious test where intersections with the faces are determined
    !-------------------------------------------------------------------------------
    Rinit(1) = SUM( Nodes % x(1:n) ) / n
    Rinit(2) = SUM( Nodes % y(1:n) ) / n
    Rinit(3) = SUM( Nodes % z(1:n) ) / n

    Hits = 0

    IF( ElemDim == 3 ) THEN
      NoFaces = BulkElement % TYPE % NumberOfFaces 
    ELSE
      NoFaces = BulkElement % TYPE % NumberOfEdges 
    END IF

    DO i=1, NoFaces          
      IF( ElemDim == 3 ) THEN
        j = BulkElement % FaceIndexes(i)
        BoundaryElement => Mesh % Faces( j )
      ELSE
        j = BulkElement % EdgeIndexes(i)
        BoundaryElement => Mesh % Edges(j)
      END IF

      CALL GetElementNodes(Nodes,BoundaryElement)

      Lambda = LineFaceIntersection2(BoundaryElement,Nodes,&
          Rinit,Rfin,Success)
      IF(.NOT. Success ) CYCLE 

      Hits = Hits + 1
      Lambdas(Hits) = Lambda
      IF (Lambda > 0.0_dp .AND. MinLambda > Lambda) MinLambda = Lambda
    END DO    

    IF( Debug ) THEN
      PRINT *,'Intersecting faces:',Hits
      PRINT *,'Lambdas:',Lambdas(1:Hits)
    END IF


    IF( Hits == 0 ) THEN
      Inside = .FALSE.
    ELSE 
      Inside = MinLambda > 1.0 - Eps .AND. MinLambda > -Eps 

      IF( Debug ) THEN
        MinDist = HUGE(MinDist)
        DO i=1,n
          Dist = SQRT( (Rfin(1)-Nodes % x(i))**2 + &
              (Rfin(2)-Nodes % y(i))**2 + &
              (Rfin(3)-Nodes % z(i))**2 )
          MinDist = MIN( Dist, MinDist )
        END DO

        PRINT *,'Dist: ', Inside,MinDist
      END IF
    END IF

  END FUNCTION SegmentElementInside


  
  !------------------------------------------------------------------------
  !> Find the particle in the mesh using actree based search. 
  !> This could be preferred in the initial finding of the correct elements.
  !> The major downside of the method is that there is no controlled face
  !> detection needed for wall interaction, for example.
  !------------------------------------------------------------------------
  SUBROUTINE LocateParticleInMeshOctree( ElementIndex, GlobalCoords, &
      LocalCoords )
    
    USE Lists
    
    INTEGER :: ElementIndex
    REAL(KIND=dp) :: GlobalCoords(3)
    REAL(KIND=dp), OPTIONAL :: LocalCoords(3)
    
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Hit, Stat
    INTEGER :: i,j,k,n
    TYPE(Nodes_t), SAVE :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Quadrant_t), POINTER, SAVE :: RootQuadrant =>NULL(), LeafQuadrant
    REAL(kind=dp) :: BoundingBox(6), eps2, eps1, uvw(3)
    
    
    Mesh => GetMesh()
    
    ! Check that the previous hit is not hit even now
    !-------------------------------------------------
    IF( ElementIndex > 0 ) THEN
      Element => Mesh % Elements( ElementIndex ) 
      n = GetElementNOFNodes(Element)
      CALL GetElementNodes(ElementNodes,Element)
      
      IF ( PointInElement( Element, ElementNodes, &
          GlobalCoords, LocalCoords ) ) RETURN
    END IF
    
    !-----------------------------------------------------------
    ! Find the right element using an octree search
    ! This is optimal when the particles are searched only once.
    !-----------------------------------------------------------
    IF ( .NOT.ASSOCIATED(Mesh % RootQuadrant) ) THEN
      BoundingBox(1) = MINVAL( Mesh % Nodes % x )
      BoundingBox(2) = MINVAL( Mesh % Nodes % y )
      BoundingBox(3) = MINVAL( Mesh % Nodes % z )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % x )
      BoundingBox(5) = MAXVAL( Mesh % Nodes % y )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
      
      eps1 = 1.0e-3
      eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
      BoundingBox(1:3) = BoundingBox(1:3) - eps2
      BoundingBox(4:6) = BoundingBox(4:6) + eps2
      
      CALL BuildQuadrantTree( Mesh,BoundingBox,Mesh % RootQuadrant)
    END IF
    RootQuadrant => Mesh % RootQuadrant
    
    Element => NULL()
    ElementIndex = 0
    CALL FindLeafElements(GlobalCoords, Mesh % MeshDim, RootQuadrant, LeafQuadrant)
    IF ( ASSOCIATED(LeafQuadrant) ) THEN
      DO i = 1, LeafQuadrant % NElemsInQuadrant
        j = LeafQuadrant % Elements(i)
        Element => Mesh % Elements(j)
        
        n = GetElementNOFNodes( Element )
        CALL GetElementNodes( ElementNodes, Element)
        
        IF ( PointInElement( Element, ElementNodes, GlobalCoords, uvw ) ) THEN
          IF( PRESENT( LocalCoords) ) LocalCoords = uvw
          ElementIndex = j
          RETURN
        END IF
      END DO
    END IF
    
    IF( ElementIndex == 0 ) THEN
      CALL Warn('LocateParticleInMeshOctree','Could not locate particle in the mesh!')
    END IF
    
  END SUBROUTINE LocateParticleInMeshOctree


  !------------------------------------------------------------------------
  !> Locate the particle using controlled marching from element to element.
  !> The crossing point between given trajectory and all face elements is 
  !> computed. The one that is passed at first is associated to the next 
  !> bulk element.
  !-------------------------------------------------------------------------
  SUBROUTINE LocateParticleInMeshMarch( ElementIndex, Rinit, Rfin, Init, &
      ParticleStatus, AccurateAtFace, StopFaceIndex, Lambda, Velo, &
      No, ParticleWallKernel, Particles )
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: ElementIndex
    REAL(KIND=dp) :: Rinit(3), Rfin(3)
    LOGICAL :: Init
    LOGICAL :: AccurateAtFace
    INTEGER :: ParticleStatus
    REAL(KIND=dp) :: Lambda
    INTEGER :: StopFaceIndex
    REAL(KIND=dp) :: Velo(3)
  
    INTEGER, OPTIONAL :: No
    OPTIONAL :: ParticleWallKernel

    INTERFACE
      SUBROUTINE ParticleWallKernel( No, r, r2, v, v2, Lambda, &
          FaceIndex, ParticleStatus ) 
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
        INTEGER :: No
        REAL(KIND=dp) :: r(3),r2(3),v(3),v2(3),Lambda
        INTEGER :: FaceIndex, ParticleStatus
      END SUBROUTINE ParticleWallKernel
    END INTERFACE

  !-------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params, BC
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Rtmp(3), Normal(3), MinLambda, eps, UnitVector(3), ds2, &
	LocalCoord(3),Velo0(3)
    LOGICAL :: Hit, DoInit, Stat, StopAtFace, AtWall, Visited = .FALSE.,&
        Debug,UseCenter,GotBC,GotBC2,ParticleBounce, Robust, Inside
    INTEGER :: i,j,k,n,dim,FaceIndex,MaxTrials,bc_id,cons_id,ElementIndex0,ParticleStatus0, &
        DebugNo, DebugPart
    INTEGER :: Problems(3), PrevNo 
    TYPE(Nodes_t), SAVE :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element, FaceElement, LeftElement, RightElement, &
        NextElement, PrevElement
    
    INTEGER, POINTER :: Neighbours(:)
    INTEGER :: NextPartition, Counter = 0
    LOGICAL, POINTER :: FaceInterface(:)
    
    SAVE :: Mesh, StopAtFace, Debug, MaxTrials, Counter, PrevNo, Eps, Robust, Problems, &
        DebugNo, DebugPart, dim
    
    Mesh => GetMesh()
    Counter = Counter + 1

    IF( .NOT. Visited ) THEN
      Params => ListGetSolverParams()
      Robust = ListGetLogical( Params,'Particle Locate Robust',stat)
      IF(.NOT. Stat) Robust = .TRUE.
      StopAtFace = ListGetLogical( Params,'Particle Stop At Face',Stat)
      MaxTrials = ListGetInteger( Params,'Max Particle Search Trials',Stat)
      IF(.NOT. Stat) MaxTrials = Mesh % NumberOfBulkElements      
     
      DebugNo = ListGetInteger( Params,'Debug particle index',Stat)
      DebugPart = ListGetInteger( Params,'Debug particle partition',Stat)

      dim = Particles % dim 
      
      Eps = ListGetConstReal( Params,'Particle Hit Tolerance',Stat)
      IF(.NOT. Stat) Eps = 1.0e-8
      Problems = 0
      Visited = .TRUE.     
    END IF

    Debug = ( No == DebugNo .AND. ParEnv % MyPe == DebugPart )

    !--------------------------------------------------------------------
    ! This is a recursive algorithm that checks the intersections 
    ! of line segments and points until correct element is found.
    ! This is optimal when the stepsize is small and there are many steps.
    !--------------------------------------------------------------------
    DoInit = Init
    IF( ElementIndex == 0 ) THEN
      DoInit = .TRUE.
      ElementIndex = 1
      UseCenter = .TRUE.
    ELSE	
      UseCenter = .NOT. AccurateAtFace
    END IF
    ElementIndex0 = ElementIndex

    IF(.NOT. UseCenter ) THEN
      ds2 = SUM( (Rinit - Rfin)**2 )
      IF( ds2 < EPSILON( ds2 ) ) RETURN
    END IF

    IF( Debug ) THEN
      PRINT *,'Starting'
      PRINT *,'Rinit:',Rinit
      PRINT *,'Rfin:',Rfin
      PRINT *,'Velo:',Velo
      PRINT *,'ds2:',ds2
    END IF
    
    NULLIFY( PrevElement ) 

    ParticleStatus = PARTICLE_LOST
    StopFaceIndex = 0
    Lambda = 1.0_dp
        
    Element => Mesh % Elements( ElementIndex ) 
    

    DO i=1,MaxTrials

      ! Use the previous element center if the true path is of no importance
      !---------------------------------------------------------------------	
      IF( UseCenter ) THEN
        n = GetElementNOFNOdes(Element)
        CALL GetElementNodes(ElementNodes,Element)
        Rtmp(1) = SUM( ElementNodes % x(1:n) ) / n
        Rtmp(2) = SUM( ElementNodes % y(1:n) ) / n
        Rtmp(3) = SUM( ElementNodes % z(1:n) ) / n
      ELSE 
        IF( i == 1 ) THEN
          Rtmp = Rinit 
        ELSE IF( .NOT. ParticleBounce ) THEN
          Rtmp = Rtmp + MinLambda * (Rfin - Rtmp) 
        END IF
      END IF

      IF( Debug ) THEN
        PRINT *,'Center:',i, UseCenter, Rtmp, MinLambda
      END IF

      
      IF( Robust ) THEN
        CALL SegmentElementIntersection2(Mesh,Element,&
            Rtmp,Rfin,MinLambda,FaceElement )
      ELSE
        CALL SegmentElementIntersection(Mesh,Element,&
            Rtmp,Rfin,MinLambda,FaceElement )
      END IF

      ParticleBounce = .FALSE.


      IF( .NOT. ASSOCIATED( FaceElement ) ) THEN
        ! One likely cause for unsuccessful operation is that the 
        ! initial node and target node are the same
        
        ds2 = SUM ( ( Rtmp - Rfin )**2 )
        IF( Debug ) THEN
          PRINT *,'NoFace:',ds2
        END IF
 
        IF( ds2 < EPSILON( ds2 ) ) THEN
          ParticleStatus = PARTICLE_HIT
          EXIT
        ELSE
!          IF( .NOT. AccurateAtFace ) THEN
!            CALL Warn('LocateParticleInMesh','No intersection found?')
!            PRINT *,'Rinit:',Rtmp
!            PRINT *,'Rfin: ',Rfin
!          END IF
          EXIT
        END IF
      ELSE IF( MinLambda > 1.0_dp - eps ) THEN
        ParticleStatus = PARTICLE_HIT
        EXIT
      ELSE 

        cons_id = FaceElement % BoundaryInfo % Constraint
        GotBC = .FALSE.
        IF ( cons_id > 0 ) THEN
          DO bc_id=1,CurrentModel % NumberOfBCs
            IF ( cons_id == CurrentModel % BCs(bc_id) % Tag ) THEN
              BC => CurrentModel % BCs(bc_id) % Values
              GotBC = .TRUE.
              EXIT
            END IF
          END DO
        END IF

        IF( Debug ) THEN
          PRINT *,'BC:',cons_id,GotBC
        END IF
 

        ! Reflect particle from face and continue the search in the same element
        !-----------------------------------------------------------------------
        GotBC2 = GotBC

        IF( GotBC ) THEN
          IF( ListGetLogical( BC,'Particle Outlet',Stat ) ) THEN
            ParticleStatus = PARTICLE_LOST
            EXIT

          ELSE IF( ListGetLogical( BC,'Particle Wall',Stat ) ) THEN
            ParticleStatus = PARTICLE_WALLBOUNDARY
            
          ELSE IF( ListGetLogical( BC,'Particle Reflect',Stat ) ) THEN
            ! First advance the particle to the point of collision
            Rtmp = Rtmp + MinLambda * (Rfin - Rtmp) 
            
            ! Then reflect the rest assuming fully elastic collision where the 
            ! normal component switches sign.
            CALL GetElementNodes(ElementNodes, FaceElement )
            Normal = NormalVector( FaceElement, ElementNodes )
            Rfin = Rfin - 2*SUM((Rfin-Rtmp)*Normal)*Normal
            
            ! Reorient the velocity vector
            UnitVector = Rfin - Rtmp
            UnitVector = UnitVector / SQRT( SUM( UnitVector** 2 ) )
            Velo = UnitVector * SQRT( SUM( Velo**2) )
                        
            IF( Debug ) THEN
              PRINT *,'Reflected',i,MinLambda,EPSILON(MinLambda)
              PRINT *,'Normal:',Normal
              PRINT *,'Rtmp:',Rtmp
              PRINT *,'Rfin:',Rfin
              PRINT *,'Abs(Velo):',SQRT(SUM(Velo**2))
              PRINT *,'Velo:',Velo
            END IF
            ParticleBounce = .TRUE.

            CYCLE

          ELSE IF( ListGetLogical( BC,'Particle Tangent',Stat ) ) THEN
            ! Get face nodes and normal vector
            CALL GetElementNodes(ElementNodes, FaceElement )
            Normal = NormalVector( FaceElement, ElementNodes, Check=.TRUE. )
            
            BLOCK
              REAL(KIND=dp) :: Amat(3,3), c(3), r0(3)

              IF( dim == 2 ) THEN
                ! Corner node of face triangle
                r0(1) = ElementNodes % x(1)
                r0(2) = ElementNodes % y(1)

                ! One basis vectors formed by the edge
                Amat(1,1) = ElementNodes % x(2) - r0(1)
                Amat(2,1) = ElementNodes % y(2) - r0(2)

                ! 2nd basis vector is the outward normal              
                Amat(1:2,2) = Normal(1:2) 

                CALL SolveLinSys2x2( Amat, c, Rfin-r0 ) 

                ! If this is outward from the normal
                IF(c(2) > 0) THEN
                  Rfin = Rfin - c(2)*Normal

                  IF( Debug ) THEN
                    PRINT *,'Tangent',i,MinLambda,EPSILON(MinLambda)
                    PRINT *,'Normal:',Normal
                    PRINT *,'Rtmp:',Rtmp
                    PRINT *,'Rfin:',Rfin
                    PRINT *,'Abs(Velo):',SQRT(SUM(Velo**2))
                    PRINT *,'Velo:',Velo
                  END IF
                  ParticleBounce = .TRUE.
                END IF

              ELSE
                ! Corner node of face triangle
                r0(1) = ElementNodes % x(1)
                r0(2) = ElementNodes % y(1)
                r0(3) = ElementNodes % z(1)
                                
                ! Two basis vectors formed by the edges
                Amat(1,1) = ElementNodes % x(2) - r0(1)
                Amat(2,1) = ElementNodes % y(2) - r0(2)
                Amat(3,1) = ElementNodes % z(2) - r0(3)

                Amat(1,2) = ElementNodes % x(3) - r0(1)
                Amat(2,2) = ElementNodes % y(3) - r0(2)
                Amat(3,2) = ElementNodes % z(3) - r0(3)

                ! 3rd basis vector is the outward normal              
                Amat(:,3) = Normal 

                CALL SolveLinSys3x3( Amat, c, Rfin-r0 ) 

                ! If this is outward from the normal
                IF(c(3) > 0) THEN
                  Rfin = Rfin - c(3)*Normal

                  IF( Debug ) THEN
                    PRINT *,'Tangent',i,MinLambda,EPSILON(MinLambda)
                    PRINT *,'Normal:',Normal
                    PRINT *,'Rtmp:',Rtmp
                    PRINT *,'Rfin:',Rfin
                    PRINT *,'Abs(Velo):',SQRT(SUM(Velo**2))
                    PRINT *,'Velo:',Velo
                  END IF
                  ParticleBounce = .TRUE.
                END IF
              END IF
            END BLOCK
            CYCLE

          ELSE IF( ListGetLogical( BC,'Particle Interact',Stat ) ) THEN
            ParticleStatus0 = ParticleStatus
            Velo0 = Velo

            CurrentModel % CurrentElement => FaceElement
            CALL ParticleWallKernel(No,Rtmp,Rfin,Velo0,Velo,MinLambda,&
                FaceElement % ElementIndex,ParticleStatus)

            ! if the particle status stays the same we are still in the same element
            IF(ParticleStatus == ParticleStatus0 ) THEN
              ParticleBounce = .TRUE.
              CYCLE
            END IF
          ELSE
            GotBC2 = .FALSE. 
          END IF
        END IF

        IF( .NOT. GotBC2 ) THEN
          LeftElement  => FaceElement % BoundaryInfo % Left
          RightElement => FaceElement % BoundaryInfo % Right
          
          
          IF( Debug ) THEN
            PRINT *,'Left and Right:',&
                ASSOCIATED(LeftElement),ASSOCIATED(RightElement)
          END IF
          
          IF( ASSOCIATED( LeftElement) .AND. ASSOCIATED(RightElement)) THEN
            IF( ASSOCIATED(Element, LeftElement)) THEN
              NextElement => RightElement
            ELSE
              NextElement => LeftElement
            END IF
            IF( StopAtFace .AND. .NOT. DoInit ) ParticleStatus = PARTICLE_FACEBOUNDARY
          ELSE
            ParticleStatus = PARTICLE_WALLBOUNDARY 
          END IF
        END IF

        ! There are different reasons why the particle is only integrated until the face
        IF( ParticleStatus == PARTICLE_WALLBOUNDARY .OR. &
            ParticleStatus == PARTICLE_FACEBOUNDARY ) THEN
                  
          Lambda = MinLambda
          StopFaceIndex = FaceElement % ElementIndex

          Rfin = Rtmp + MinLambda * (Rfin - Rtmp) 
          Velo = 0.0_dp

          IF( Debug ) THEN
            PRINT *,'WallBC:',Rfin, MinLambda
          END IF

          EXIT                      
        END IF
      END IF

      IF( Debug ) THEN
        PRINT *,'Same Elements:', ASSOCIATED( NextElement ), &
            ASSOCIATED( NextElement, Element), ASSOCIATED( NextElement, PrevElement )
      END IF
      
      ! continue the search to new elements
      IF( .NOT. ASSOCIATED( NextElement ) ) THEN
        CALL Warn('LocateParticleInMeshMarch','Element not associated!')
      END IF

      IF( ASSOCIATED( NextElement, Element ) ) THEN
        CALL Warn('LocateParticleInMeshMarch','Elements are the same!')
      END IF


      IF( ASSOCIATED( NextElement, PrevElement ) ) THEN
        CALL GetElementNodes(ElementNodes,NextElement)
        IF( Robust ) THEN
          Inside =  SegmentElementInside(Mesh,NextElement,Rfin,Debug) 
        ELSE
          Inside = PointInElement( NextElement, ElementNodes, Rfin, LocalCoord )
        END IF

        IF( Debug ) THEN
          PRINT *,'NextElement',Inside,Robust,NextElement % ElementIndex,LocalCoord
        END IF

        IF( Inside ) THEN
	  Problems(1) = Problems(1) + 1
           CALL Warn('LocateParticleInMeshMarch','Elements are same, found in NextElement!')          
          Element => NextElement
          ParticleStatus = PARTICLE_HIT	  
          EXIT
        END IF 

        CALL GetElementNodes(ElementNodes,Element)
        IF( Robust ) THEN
          Inside =  SegmentElementInside(Mesh,Element,Rfin,Debug) 
        ELSE
          Inside = PointInElement( Element, ElementNodes, Rfin, LocalCoord )
        END IF

        IF( Debug ) THEN
          PRINT *,'ThisElement',Inside,Robust,Element % ElementIndex,LocalCoord
        END IF

       IF( Inside ) THEN
	  Problems(2) = Problems(2) + 1
           CALL Warn('LocateParticleInMeshMarch','Elements are same, found in Element!')          
          ParticleStatus = PARTICLE_HIT	  
          EXIT
        END IF 

        Problems(3) = Problems(3) + 1
        WRITE(Message,'(A,3ES12.3)') 'Losing particle '//I2S(No)//' in: ',Rfin(1:3)
        CALL Info('LocateParticlesInMesh',Message,Level=15)
        
        ParticleStatus = PARTICLE_LOST
        EXIT
      END IF

      PrevElement => Element
      Element => NextElement
    END DO

    IF( i >= MaxTrials ) THEN
      PRINT *,'Used maximum number of trials',MaxTrials,No,SQRT(ds2)
    END IF

    IF( ParticleStatus == PARTICLE_LOST ) THEN
      ElementIndex = 0
    ELSE	
      ElementIndex = Element % ElementIndex
    END IF

100 CONTINUE

    ! This is just for debugging 
    IF( No < PrevNo .AND. Problems(3) > 0 ) THEN
      WRITE( Message,'(A,3I0)') 'Problems in locating particles:',Problems
      CALL Info('LocateParticleInMeshMarch',Message,Level=10)
    END IF
    Problems = 0
    PrevNo = No
     
  END SUBROUTINE LocateParticleInMeshMarch



  !------------------------------------------------------------------------
  !> Locate the particles in their new positions in the mesh.
  !-------------------------------------------------------------------------
  SUBROUTINE LocateParticles( Particles, ParticleWallKernel )
    
    USE Lists
    
    TYPE(Particle_t), POINTER :: Particles
    OPTIONAL :: ParticleWallKernel

    INTERFACE 
      SUBROUTINE ParticleWallKernel( No, r, r2, v, v2, Lambda, &
          FaceIndex, ParticleStatus ) 
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
        INTEGER :: No
        REAL(KIND=dp) :: r(3),r2(3),v(3),v2(3),Lambda
        INTEGER :: FaceIndex, ParticleStatus
      END SUBROUTINE ParticleWallKernel
    END INTERFACE


    !-----------------------------------------------------------------------
    LOGICAL :: PartitionChangesOnly 
    INTEGER :: PartitionChanges, Status, ElementIndex, No, &
               NoParticles, dim, ElementIndex0
    REAL(KIND=dp) :: Rinit(3), Rfin(3),Rfin0(3),Velo(3), Velo0(3), dtime, Speed
    LOGICAL :: Stat, InitLocation, AccurateAtFace, AccurateAlways, AccurateNow, debug
    INTEGER :: FaceIndex, FaceIndex0, Status0, InitStatus, MaxIter, Iter
    REAL(KIND=dp) :: Lambda
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: DtVar
    CHARACTER(*), PARAMETER :: Caller = 'LocateParticles'

    CALL Info(Caller,'Locating particles in mesh',Level=10)

    Params => ListGetSolverParams()
    Mesh => GetMesh()
    dim = Particles % dim      
    PartitionChangesOnly = .FALSE.
    Velo = 0.0_dp
    Debug = .FALSE.

    ! Particles may be located either using the mid-point of the current element as the 
    ! reference point for finding face intersections, or using the true starting 
    ! point which is more prone to epsilon-errors.
    !---------------------------------------------------------------------------
    AccurateAlways = ListGetLogical( Params,'Particle Accurate Always',Stat)
    AccurateAtFace = ListGetLogical( Params,'Particle Accurate At Face',Stat)

    
    IF( .NOT. Particles % DtConstant ) THEN
      DtVar => ParticleVariableGet( Particles,'particle dt')
      IF(.NOT. ASSOCIATED( DtVar ) ) THEN
        CALL Fatal(Caller,'Variable timestep, > particle dt < should exist!')
      END IF
    END IF

    MaxIter = 10
    Iter = 0

100 NoParticles = Particles % NumberOfParticles
    Iter = Iter + 1

    CALL Info(Caller,'Locating particles iteration: '//I2S(Iter),Level=12)
    
    !Debug = ( iter > 8 ) 
    
    DO No = 1, NoParticles
      
      Status = Particles % Status( No )
      InitStatus = Status
   
      IF( Status >= PARTICLE_LOST ) CYCLE
      IF( Status < PARTICLE_INITIATED ) CYCLE
      IF( Status == PARTICLE_WALLBOUNDARY ) CYCLE
      IF( Status == PARTICLE_FIXEDCOORD ) CYCLE

      IF( .NOT. Particles % DtConstant ) THEN
        IF( ABS( DtVar % Values(No) ) < TINY( dtime ) ) CYCLE
      END IF
              
      IF ( PartitionChangesOnly .AND. Status /= PARTICLE_PARTBOUNDARY ) CYCLE
      
      InitLocation = ( Status < PARTICLE_LOCATED ) 
      Velo = GetParticleVelo( Particles, No )       
      IF( Status == PARTICLE_INITIATED ) THEN
        AccurateNow = .FALSE.
      ELSE
        AccurateNow = AccurateAlways
      END IF
      FaceIndex0 = 0
      ElementIndex0 = 0
      
200   ElementIndex = GetParticleElement( Particles, No )
      Rfin = GetParticleCoord( Particles, No )
      Velo = GetParticleVelo( Particles, No )      
      IF( AccurateNow ) Rinit = GetParticlePrevCoord( Particles, No )        
      Rinit = GetParticlePrevCoord( Particles, No )        
      Speed = SQRT(SUM(Velo**2))
      
      IF( debug ) THEN
        PRINT *,parenv % mype, 'going No',No,'Element',ElementIndex,'Face',FaceIndex,'Status',Status
        PRINT *,parenv % mype, 'going Init:    ',Rinit(1:dim),Rfin(1:dim)
        PRINT *,parenv % mype, 'going Velo:',GetParticleVelo(Particles,No), Velo(1:dim)
        PRINT *,parenv % mype, 'goint Speed:',speed, Status
      END IF

      ! If the particle speed is located and its speed is zero this is not going anywhere...
      IF(Speed < EPSILON(Speed) ) THEN
        IF( Status > PARTICLE_INITIATED ) CYCLE
      END IF
        
      CALL LocateParticleInMeshMarch(ElementIndex, Rinit, Rfin, InitLocation, &
          Status,AccurateNow, FaceIndex, Lambda, Velo, No, ParticleWallKernel, Particles )

      IF( debug ) THEN
        PRINT *,parenv % mype, 'leaving No',No,'Element',ElementIndex,'Face',FaceIndex,'Status',Status
        PRINT *,parenv % mype, 'leaving Init:    ',Rinit(1:dim),Rfin(1:dim)
        PRINT *,parenv % mype, 'Velo:',GetParticleVelo(Particles,No), Velo(1:dim)
      END IF

      ! If a boundary face is passed then repeat the process more diligently. 
      ! Note that if we want proper collisions they should be implemented below
      !----------------------------------------------------------------------------
      IF( .NOT. AccurateNow ) THEN
        AccurateNow = AccurateAtFace .AND. FaceIndex > 0

        IF ( debug ) PRINT*,accuratenow, accurateatface, faceindex > 0

        IF( AccurateNow ) THEN
          FaceIndex0 = FaceIndex
          ElementIndex0 = ElementIndex
          Status0 = Status
          Rfin0 = Rfin        
          Velo0 = Velo
          ElementIndex = GetParticleElement( Particles, No )
          Rfin = GetParticleCoord( Particles, No )
          Velo = GetParticleVelo( Particles, No )      
          IF ( debug ) PRINT*,parenv % mype, 'go 200 '; FLUSH(6)
          GOTO 200
        END IF
      END IF

      IF( FaceIndex0 > 0 .AND. FaceIndex == 0 ) THEN        
        ! Currently it is assumed that if success with the second method is not obtained then
        ! the particle is already sitting on the face observed by the more robust method. 
        !-------------------------------------------------------------------------------------
        IF( .FALSE. ) THEN
          CALL Warn(Caller,'Difference between robust and accurate?')
        END IF   

        Status = Status0
        Rfin = Rfin0
        Velo = Velo0
        ElementIndex = ElementIndex0
        FaceIndex = FaceIndex0
      END IF


      IF( ElementIndex == 0 ) THEN
        Rfin = Rfin0
        Velo = Velo0
        ElementIndex = ElementIndex0
        FaceIndex = FaceIndex0
      END IF
      
      IF( debug ) THEN
        PRINT *,parenv % mype, 'No',No,'Element',ElementIndex,'Face',FaceIndex,'Status',Status
        PRINT *,parenv % mype, 'Init:    ',Rinit(1:dim),Rfin(1:dim)
        PRINT *,parenv % mype, 'Velo:',GetParticleVelo(Particles,No), Velo(1:dim)
      END IF

      Particles % ElementIndex(No) = ElementIndex       
      Particles % Status(No) = Status
      Particles % FaceIndex(No) = FaceIndex

      CALL SetParticleCoord( Particles, No, Rfin  )                
      IF( ElementIndex == 0 ) Velo = 0.0_dp

      CALL SetParticleVelo( Particles, No, Velo  )
    END DO

    ! Change the partition in where the particles are located
    ! Only applies to parallel cases.
    !------------------------------------------------------------------------
    PartitionChanges = ChangeParticlePartition( Particles )
    
    IF( PartitionChanges > 0 .AND. Iter < MaxIter ) THEN      
      CALL Info(Caller,'Parallel locate loop '//I2S(iter)//' with '&
          //I2S(PartitionChanges)//' particles!',Level=7)
      PartitionChangesOnly = .TRUE.
      GOTO 100
    END IF
    
  END SUBROUTINE LocateParticles
  
  
  
  !--------------------------------------------------------------------------
  !> Given the element & global coordinates returns the local coordinates.
  !> The idea of this routine is to transparently block the local coordinate
  !> search from the user by directly giving the basis function values related
  !> to a global coordinate. Sloppy tolerances are used since we *should* 
  !> have already located the element.
  !--------------------------------------------------------------------------
  FUNCTION ParticleElementInfo( Element, GlobalCoord, &
      SqrtElementMetric, Basis, dBasisdx ) RESULT ( stat )
    
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: GlobalCoord(:), SqrtElementMetric, LocalDistance
    REAL(KIND=dp) :: Basis(:)
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)
    LOGICAL :: Stat, Debug
    INTEGER :: Misses(2) = 0
    CHARACTER(*), PARAMETER :: Caller = 'ParticleElementInfo'
    
    SAVE Misses    

    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: LocalCoord(3),u,v,w
    INTEGER :: n
    
    SAVE ElementNodes

    IF(.NOT. ASSOCIATED(Element) ) THEN
      CALL Fatal('ParticleElementInfo','Element not associated!')
    END IF

    IF(.NOT. ASSOCIATED(Element % TYPE) ) THEN
      PRINT *,'Element:',Element % ElementIndex
      CALL Fatal('ParticleElementInfo','Element % Type not associated!')
    END IF
    
    n = Element % TYPE % NumberOfNodes
    CALL GetElementNodes(ElementNodes,Element)
    
    Stat = PointInElement( Element, ElementNodes, &
        GlobalCoord, LocalCoord, GlobalEps = -1.0_dp, LocalEps = 1.0e3_dp, &
	LocalDistance = LocalDistance ) 

    IF( .NOT. Stat ) THEN
      Misses(1) = Misses(1) + 1

      IF( MODULO( SUM( Misses ), 101 ) == 100 ) PRINT *,'Misses:',Misses

      IF( .FALSE.) THEN
        IF( .NOT. Stat ) THEN
          CALL Warn(Caller,'Should have found the node!')
        ELSE
          CALL Warn(Caller,'Distance from element higher than expected!')
        END IF
        PRINT *,'LocalDistance:',LocalDistance,'Element:',Element % ElementIndex
        PRINT *,'Nodes X:',ElementNodes % x(1:n) - GlobalCoord(1)
        PRINT *,'Nodes Y:',ElementNodes % y(1:n) - GlobalCoord(2)
        PRINT *,'Nodes Z:',ElementNodes % z(1:n) - GlobalCoord(3)
      END IF
      RETURN
    END IF

    u = LocalCoord(1)
    v = LocalCoord(2)
    w = LocalCoord(3)
    
    stat = ElementInfo( Element, ElementNodes, U, V, W, SqrtElementMetric, &
        Basis, dBasisdx )
    IF(.NOT. Stat) Misses(2) = Misses(2) + 1
    
  END FUNCTION ParticleElementInfo
  


  !-------------------------------------------------------------------------
  !> The routine returns velocity and optionally a gradient of velocity.
  !> These kind of functions are needed repeated and therefore to reduced the 
  !> size of individual solvers it has been hard coded here. 
  !--------------------------------------------------------------------------
  
  SUBROUTINE GetVectorFieldInMesh(Var, Element, Basis, Velo, dBasisdx, GradVelo )
    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t) :: Element
    REAL(KIND=dp) :: Basis(:), Velo(:) 
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:), GradVelo(:,:)
    
    TYPE(Valuelist_t), POINTER :: Params
    INTEGER, POINTER :: LocalPerm(:)
    REAL(KIND=dp), POINTER :: LocalVelo(:,:)
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: Dofs
    REAL(KIND=dp) :: SumBasis
    INTEGER :: i,j,k,n,npos,ind,dim
    LOGICAL :: GotIt, InterfaceNodes
    LOGICAL :: Visited
    
    
    SAVE :: Visited, Dim, LocalVelo, LocalPerm, InterfaceNodes
    
    IF(.NOT. Visited ) THEN
      Mesh => GetMesh()
      Params => ListGetSolverParams()
      n = Mesh % MaxElementNodes
      ALLOCATE( LocalPerm(n), LocalVelo(n,3) )
      
      InterfaceNodes = GetLogical( Params,'Interface Nodes',GotIt)

      LocalPerm = 0
      LocalVelo = 0.0_dp
      Dim = Mesh % MeshDim
      Visited = .TRUE.
    END IF

    Velo = 0.0_dp
    IF( PRESENT( GradVelo ) ) GradVelo = 0.0_dp

    IF(.NOT. ASSOCIATED(Var) ) THEN
      CALL Fatal('GetVectorFieldInMesh','Variable not associated!')
    END IF

    
    n = Element % TYPE % NumberOfNodes    
    IF( Var % TYPE == Variable_on_nodes_on_elements ) THEN      
      LocalPerm(1:n) = Var % Perm( Element % DGIndexes )
    ELSE
      LocalPerm(1:n) = Var % Perm( Element % NodeIndexes )
    END IF
    npos = COUNT ( LocalPerm(1:n) > 0 )
        
    
    IF( npos == 0 ) RETURN

    ! MAX 3 velocity direction (4th variable is pressure)
    ! This has to be accounted for in the Var % Values permutations
    dofs = MIN(3,Var % Dofs)
    
    !-----------------------------------------------------------------
    ! compute the velocity also for case when the particle
    ! has just crossed the boundary. For example, its floating on the 
    ! fluid boundary. This is a little bit fishy and could perhaps 
    ! only be done conditionally....
    ! Can't really determine the gradient here
    !-----------------------------------------------------------------
    IF( npos == n ) THEN
      DO i=1,n
        j = LocalPerm(i)
	DO k=1,dofs
          ! For correct permutation we have to account for the full dimension
          ! of the Flow Solution (e.g., Stokes 3D : vx, vy,vz, p)
          LocalVelo(i,k) = Var % Values( Var % Dofs * (j-1) + k)
        END DO
      END DO
    ELSE    
      IF(.NOT. InterfaceNodes ) RETURN

      SumBasis = 0.0_dp
      DO i=1,n
        j = LocalPerm(i)
        IF( j > 0 ) THEN
          SumBasis = SumBasis + Basis(i)
          DO k=1,dofs
            LocalVelo(i,k) = Var % Values( Var % Dofs * (j-1) + k)
          END DO
        ELSE
          Basis(i) = 0.0_dp
          LocalVelo(i,1:dim) = 0.0_dp
        END IF
      END DO
    END IF
    

    DO i=1,dofs
      Velo(i) = SUM( Basis(1:n) * LocalVelo(1:n,i) )
      IF( PRESENT( GradVelo ) ) THEN
        ! dBasisdx has only stuff up until the dimension!
        DO j=1,dim
          GradVelo(i,j) = SUM( dBasisdx(1:n,j) * LocalVelo(1:n,i) )
        END DO
      END IF
    END DO
    
    IF( npos < n ) THEN
      Velo(1:dofs) = Velo(1:dofs) / SumBasis
      IF( PRESENT( GradVelo ) ) THEN
        GradVelo(:,1:dim) = GradVelo(:,1:dim) / SumBasis
      END IF
    END IF
   
  END SUBROUTINE GetVectorFieldInMesh


  !-------------------------------------------------------------------------
  !> The routine returns a potential and its gradient.
  !--------------------------------------------------------------------------
  
  SUBROUTINE GetScalarFieldInMesh(Var, Element, Basis, Pot, dBasisdx, GradPot )
    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t) :: Element
    REAL(KIND=dp) :: Basis(:), Pot 
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:), GradPot(:)
    
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: LocalPerm(:)
    REAL(KIND=dp), POINTER :: LocalField(:)
    INTEGER :: i,j,n,dim
    LOGICAL :: Visited
    
    
    SAVE :: Visited, Mesh, Dim, LocalPerm, LocalField
    
    IF(.NOT. Visited ) THEN
      Mesh => GetMesh()
      n = Mesh % MaxElementNodes
      ALLOCATE( LocalPerm(n), LocalField(n) )
      LocalPerm = 0
      LocalField = 0.0_dp
      Dim = Mesh % MeshDim
      Visited = .TRUE.
    END IF
     
    Pot = 0.0_dp
    IF( PRESENT( GradPot ) ) GradPot = 0.0_dp
    
    IF(.NOT. ASSOCIATED( Var ) ) RETURN
    
    n = Element % TYPE % NumberOfNodes
    IF( ASSOCIATED( Var % Perm ) ) THEN
      LocalPerm(1:n) = Var % Perm( Element % NodeIndexes )
      IF( .NOT. ALL ( LocalPerm(1:n) > 0 )) RETURN
      LocalField(1:n) = Var % Values( LocalPerm(1:n) )
    ELSE
      ! Some variables do not have permutation, most importantly the node coordinates
      LocalField(1:n) = Var % Values( Element % NodeIndexes )
    END IF

    Pot = SUM( Basis(1:n) * LocalField(1:n) )

    IF( PRESENT( GradPot ) ) THEN
      DO i=1,dim
        GradPot(i) = SUM( dBasisdx(1:n,i) * LocalField(1:n) )
      END DO
    END IF
    
  END SUBROUTINE GetScalarFieldInMesh


  
  !-------------------------------------------------------------------------
  !> The routine returns the possible intersection of a secondary element 
  !> with a different material property and the circle / sphere.
  !> For example, the buoyancy at the interface will depend on the weighted
  !> sum of the densities of the two materials. 
  !--------------------------------------------------------------------------
  
  FUNCTION GetParticleElementIntersection(Particles,BulkElement, Basis, Coord, &
      Radius, BulkElement2, VolumeFraction, AreaFraction ) RESULT ( Intersect )
    
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Element_t), POINTER :: BulkElement, BulkElement2
    REAL(KIND=dp) :: Basis(:)
    REAL(KIND=dp) :: Coord(3), Radius, VolumeFraction
    REAL(KIND=dp), OPTIONAL :: AreaFraction
    LOGICAL :: Intersect
    
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Dist, Normal(3), SumBasis
    TYPE(ValueList_t), POINTER :: Material, Material2, BC
    TYPE(Element_t), POINTER :: BoundaryElement, Left, Right
    TYPE(Nodes_t) :: BoundaryNodes
    INTEGER :: i,j,k,n,imax,body_id,body_id2,mat_id,mat_id2,dim,ind
    LOGICAL :: Visited
    
    
    SAVE :: Visited, Mesh, Dim
    
    IF(.NOT. Visited ) THEN
      Mesh => GetMesh()
      Dim = Mesh % MeshDim
      Visited = .TRUE.    
    END IF
    
    Intersect = .FALSE.
    VolumeFraction = 0.0_dp
    
    ! This element has no boundary / material interface
    IF( Particles % InternalElements( BulkElement % ElementIndex ) ) RETURN
    
    ! If the radius of the particle is zero then it sees only the properties of one point
    IF( Radius < TINY( Radius ) ) RETURN
    
    n = BulkElement % TYPE % NumberOfNodes
    body_id = BulkElement % BodyId
    mat_id = ListGetInteger( CurrentModel % Bodies(body_id) % Values,'Material' )
    
    IF( dim == 3 ) THEN
      imax = BulkElement % TYPE % NumberOfFaces 
    ELSE
      imax = BulkElement % TYPE % NumberOfEdges  
    END IF
    
    DO i=1, imax
      
      IF( dim == 3 ) THEN
        j = BulkElement % FaceIndexes(i)
        BoundaryElement => Mesh % Faces( j )
      ELSE
        j = BulkElement % EdgeIndexes(i)
        BoundaryElement => Mesh % Edges(j)
      END IF
      
      IF( .NOT. ASSOCIATED( BoundaryElement % BoundaryInfo ) ) CYCLE
      
      Left => BoundaryElement % BoundaryInfo % Left
      Right => BoundaryElement % BoundaryInfo % Right
      
      IF(.NOT. (ASSOCIATED( Left ) .AND. ASSOCIATED( Right ) ) ) CYCLE
      
      IF( ASSOCIATED( BulkElement, Right ) ) THEN
        BulkElement2 => Left
      ELSE 
        BulkElement2 => Right 
      END IF
      
      IF( .NOT. ASSOCIATED( BulkElement2 ) ) CYCLE
      
      body_id2 = BulkElement2 % BodyId
      
      IF( body_id2 > CurrentModel % NumberOfBodies ) THEN
        PRINT *,'BodyIds:',body_id,body_id2,CurrentModel % NumberOfBodies
        PRINT *,'ElemIds:',BulkElement % ElementIndex, BulkElement2 % ElementIndex
        PRINT *,'Types:',BulkElement % TYPE % NumberOfNodes, &
            BulkElement2 % TYPE % NumberOfNodes
        body_id2 = 0
      END IF
      
      IF( body_id2 == 0 ) CYCLE
      
      mat_id2 = ListGetInteger( CurrentModel % Bodies(body_id2) % Values,'Material' )
      
      ! If the materials are the same the density is ok
      IF( mat_id2 == mat_id ) CYCLE          

      ! If there is an material interface, check for distance
      CALL GetElementNodes(BoundaryNodes,BoundaryElement)
      Dist = PointFaceDistance(BoundaryElement,BoundaryNodes,Coord,Normal)
      Dist = ABS( Dist )       
      
      ! Is is assumed that each element may only have one density interface
      IF( Dist > Radius ) RETURN
      
      IF( dim == 3 ) THEN
        ! based on the formula of sphere-sphere intersection as in Wolfram MathWorld
        VolumeFraction = (Radius + Dist / 2 ) * (Radius - Dist)**2 / Radius**3
        IF( PRESENT( AreaFraction ) ) THEN
          AreaFraction = ( 1.0_dp - Dist/Radius )/2.0_dp
        END IF
      ELSE
        ! based on the formula of circle-circle intersection as in Wolfram MathWorld
        VolumeFraction = ( ( Radius ** 2) * ACOS( Dist / Radius ) &
            - Dist * SQRT( Radius ** 2 - Dist ** 2 ) ) / (PI * Radius**2) 
        IF( PRESENT( AreaFraction ) ) THEN
          AreaFraction = ACOS( Dist / Radius ) / PI
        END IF
      END IF
      
      !     PRINT *,'VolumeFraction:',VolumeFraction, Density
      RETURN
    END DO
    
  END FUNCTION GetParticleElementIntersection


  !-------------------------------------------------------------
  !> This subroutine may be used to enquire position dependent material data.
  !> Also if the particle is split between two elements then this 
  !> routine can assess the data on the secondary mesh.
  !-------------------------------------------------------------
  FUNCTION GetMaterialPropertyInMesh(PropertyName, BulkElement, Basis, &
      BulkElement2, VolumeFraction ) RESULT ( Property )
    
    CHARACTER(LEN=*) :: PropertyName
    TYPE(Element_t), POINTER :: BulkElement
    REAL(KIND=dp) :: Basis(:)
    TYPE(Element_t), POINTER, OPTIONAL :: BulkElement2
    REAL(KIND=dp), OPTIONAL :: VolumeFraction
    REAL(KIND=dp) :: Property
    
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp), POINTER :: ElemProperty(:)
    REAL(KIND=dp) :: Property2
    TYPE(ValueList_t), POINTER :: Material, Material2
    INTEGER :: i,j,k,n,mat_id,mat_id2
    LOGICAL :: Visited
    
    
    SAVE :: Visited, Mesh, ElemProperty
    
    IF(.NOT. Visited ) THEN
      Mesh => GetMesh()
      n = Mesh % MaxElementNodes
      ALLOCATE( ElemProperty( n ) )
      ElemProperty = 0.0_dp
      Visited = .TRUE.    
    END IF
    
    NodeIndexes => BulkElement % NodeIndexes
    n = BulkElement % TYPE % NumberOfNodes
    mat_id = ListGetInteger( CurrentModel % Bodies(BulkElement % BodyId) % Values,'Material' )
    Material => CurrentModel % Materials(mat_id) % Values
    
    ElemProperty(1:n) = ListGetReal( Material,PropertyName,n,NodeIndexes) 
    Property = SUM( Basis(1:n) * ElemProperty(1:n) )
    
    IF( .NOT. PRESENT ( VolumeFraction ) ) RETURN
    IF( .NOT. PRESENT ( BulkElement2 ) ) RETURN
    IF( VolumeFraction < TINY( VolumeFraction) ) RETURN
    
    IF( ASSOCIATED( BulkElement2 ) ) THEN
      mat_id2 = ListGetInteger( CurrentModel % Bodies(BulkElement2 % BodyId) % Values,'Material' )
    ELSE
      mat_id2 = 0
    END IF
    
    ! If the materials are the same the density is ok
    IF( mat_id2 == mat_id ) RETURN
    
    ! If there is an material interface, check for distance
    IF( mat_id2 == 0 ) THEN
      Property2 = 0.0_dp
    ELSE
      NodeIndexes => BulkElement2 % NodeIndexes
      n = BulkElement2 % TYPE % NumberOfNodes
      Material2 => CurrentModel % Materials(mat_id2) % Values
      
      ElemProperty(1:n) = ListGetReal( Material,PropertyName,n,NodeIndexes) 
      
      ! One cannot use the basis functions of the primary element. 
      ! and this is valid for cases with constant material parameters.
      !------------------------------------------------------------------
      Property2 = SUM( ElemProperty(1:n) ) / n
    END IF
    
    Property = VolumeFraction * Property2 + (1-VolumeFraction) * Property
    
  END FUNCTION GetMaterialPropertyInMesh


  !-------------------------------------------------------------
  !> This routine creates the nearest neighbours for all nodes.
  !> The particle-particle connections may then be found by going
  !> through all the nodes of elements.
  !-------------------------------------------------------------
  SUBROUTINE CreateNeighbourList( Particles ) 
    
    TYPE(Particle_t), POINTER :: Particles
    
    INTEGER :: ElementIndex, dim
    REAL(KIND=dp) :: Coord(3), dist, mindist
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: i,j,k,n,node
    TYPE(Nodes_t), SAVE :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: NoNodes, NoParticles, MaxClosest  
    
    Mesh => GetMesh()
    NoNodes = Mesh % NumberOfNodes
    NoParticles = Particles % NumberOfParticles
    dim = Particles % dim 
    
    IF( .NOT. Particles % NeighbourTable ) THEN
      ALLOCATE( Particles % NoClosestParticle( NoNodes ) ) 
      ALLOCATE( Particles % CumClosestParticle( NoNodes+1 ) ) 
      ALLOCATE( Particles % ClosestNode(NoParticles) )
      Particles % NeighbourTable = .TRUE.
    END IF
    
    IF( SIZE( Particles % ClosestNode ) < NoParticles ) THEN
      CALL Fatal('CreateNeighbourList','ClosestNode vector of wrong size')
    END IF
    
    ! First find the closest node to each particle
    !-----------------------------------------------
    Particles % ClosestNode = 0
    Particles % NoClosestParticle = 0
    DO i=1,NoParticles 
      IF( Particles % Status(i) >= PARTICLE_LOST ) CYCLE
      IF( Particles % Status(i) < PARTICLE_INITIATED ) CYCLE
      
      ElementIndex = Particles % ElementIndex(i)
      Element => Mesh % Elements( ElementIndex )
      n = GetElementNOFNodes(Element)
      CALL GetElementNodes(ElementNodes,Element)    
      Coord(1:dim) = Particles % Coordinate(i,1:dim)
      
      ! Find the minimum distance node (using squares is faster)
      mindist = HUGE( mindist ) 
      DO j=1,n
        dist = ( ElementNodes % x(j) - Coord(1) )**2
        dist = dist +  ( ElementNodes % y(j) - Coord(2) )**2
        IF( dim == 3 ) THEN
          dist = dist +  ( ElementNodes % z(j) - Coord(3) )**2
        END IF
        IF( dist < mindist ) THEN
          mindist = dist 
          k = j
        END IF
      END DO
      node = Element % NodeIndexes(k)
      Particles % ClosestNode(i) = node
      Particles % NoClosestParticle(node) = Particles % NoClosestParticle(node) + 1
    END DO


    ! For parallel computation create a secondary copy of the neighbouring particles
    ! marked as ghost particle and update the total number.
    !-------------------------------------------------------------------------------
    CALL CreateGhostParticles( Particles )

    Particles % FirstGhost = NoParticles + 1
    IF( Particles % NumberOfParticles > NoParticles ) THEN
      Particles % FirstGhost = NoParticles + 1
      NoParticles = Particles % NumberOfParticles
    END IF
    
    ! Count the cumulative number of closest particles for given node
    !-----------------------------------------------------------------
    Particles % CumClosestParticle(1) = 1
    MaxClosest = 0
    DO i=1,NoNodes
      j = Particles % NoClosestParticle(i)
      MaxClosest = MAX( MaxClosest, j )
      Particles % CumClosestParticle(i+1) = Particles % CumClosestParticle(i)+j
    END DO
    Particles % MaxClosestParticles = MaxClosest

    ! And finally, add the closest neighbors to the table
    !----------------------------------------------------------------
    IF ( ASSOCIATED(Particles % ClosestParticle) ) &
        DEALLOCATE(Particles % ClosestParticle )
    ALLOCATE( Particles % ClosestParticle(Particles % CumClosestParticle(NoNodes+1)) )
    
    Particles % NoClosestParticle = 0
    Particles % ClosestParticle = 0
    DO i=1,NoParticles     
      IF ( Particles % Status(i) == PARTICLE_LOST ) CYCLE
      IF ( Particles % Status(i) < PARTICLE_INITIATED ) CYCLE
      node = Particles % ClosestNode(i) 
      j = Particles % NoClosestParticle(node) 
      k = Particles % CumClosestParticle(node)
      Particles % ClosestParticle(k+j) = i
      Particles % NoClosestParticle(node) = j + 1
    END DO
    
  END SUBROUTINE CreateNeighbourList
  

  SUBROUTINE CreateGhostParticles(Particles)
    TYPE(Particle_t), POINTER :: Particles
    !---------------------------------------------------------
    INTEGER i,j,k,l,m,n,dim,NoPartitions, node, &
        Proc, ierr, status(MPI_STATUS_SIZE), n_part, nReceived
    
    INTEGER, ALLOCATABLE :: Perm(:), Indexes(:), Neigh(:), &
        Recv_parts(:), Requests(:)
    TYPE(Mesh_t), POINTER :: Mesh
    
    TYPE(ParallelInfo_t), POINTER :: PI
    
    LOGICAL, ALLOCATABLE :: IsNeighbour(:)
    INTEGER, POINTER :: Neighbours(:), Closest(:)
    
    TYPE ExchgInfo_t
      INTEGER :: n=0
      INTEGER, ALLOCATABLE :: Gindex(:), Particles(:)
    END TYPE ExchgInfo_t
    
    REAL(KIND=dp), ALLOCATABLE :: Buf(:)
    TYPE(ExchgInfo_t), POINTER :: Info(:)
    !--------------------------------------------------------
    
    nReceived = 0
    IF( ParEnv% PEs == 1 ) RETURN
    
    Mesh => GetMesh()
    dim = Particles % dim
    
    ! Count & Identify neighbouring partitions:
    ! -----------------------------------------
    ALLOCATE(IsNeighbour(ParEnv % PEs))
    NoPartitions = MeshNeighbours(Mesh,IsNeighbour)
    ALLOCATE(Perm(ParEnv % PEs), Neigh(NoPartitions) )
    Perm = 0
    
    NoPartitions=0
    DO i=1,ParEnv % PEs
      IF ( i-1==ParEnv % Mype ) CYCLE
      IF ( IsNeighbour(i) ) THEN
        NoPartitions=NoPartitions+1
        Perm(i) = NoPartitions
        Neigh(NoPartitions) = i-1
      END IF
    END DO
    DEALLOCATE(IsNeighbour)
    
    ! Receive interface sizes:
    !--------------------------
    ALLOCATE( Recv_Parts(NoPartitions), Requests(NoPartitions) )
    DO i=1,NoPartitions
      CALL MPI_iRECV( Recv_Parts(i),1, MPI_INTEGER, Neigh(i), &
          2000, ELMER_COMM_WORLD, requests(i), ierr )
    END DO
    
    PI => Mesh % ParallelInfo
    
    ! Exchange interface particles
    ! ----------------------------
    ALLOCATE(Info(NoPartitions))
    DO i=1,NoPartitions
      Info(i) % n = 0
    END DO
    
    DO i=1,Particles % NumberOfParticles
      IF ( Particles % Status(i) == PARTICLE_LOST ) CYCLE
      
      node = Particles % ClosestNode(i)
      IF ( .NOT. PI % GInterface(node) ) CYCLE
      Neighbours => PI % NeighbourList(node) % Neighbours
      DO j=1,SIZE(Neighbours)
        proc = Neighbours(j)
        IF ( Proc==Parenv % mype ) CYCLE
        proc = Perm(proc+1)
        IF ( Proc<=0 ) CYCLE
        Info(proc) % n = Info(proc) % n+1
      END DO
    END DO
    
    DO i=1,NoPartitions
      CALL MPI_BSEND( Info(i) % n, 1, MPI_INTEGER, Neigh(i), &
          2000, ELMER_COMM_WORLD, ierr )
    END DO
    
    !
    ! Collect particles to be sent to neighbours:
    ! -------------------------------------------
    DO i=1,NoPartitions
      IF ( Info(i) % n==0 ) CYCLE
      ALLOCATE( Info(i) % Gindex(Info(i) % n), Info(i) % Particles(Info(i) % n) )
      Info(i) % n = 0
    END DO
    
    DO i=1,Particles % NumberOfParticles
      IF ( Particles % Status(i) == PARTICLE_LOST ) CYCLE
      
      node = Particles % ClosestNode(i)
      IF ( .NOT. PI % GInterface(node) ) CYCLE
      Neighbours => PI % NeighbourList(node) % Neighbours
      DO j=1,SIZE(Neighbours)
        proc = Neighbours(j)
        IF ( Proc==Parenv % mype ) CYCLE
        proc = Perm(proc+1)
        IF ( Proc<=0 ) CYCLE
        Info(proc) % n = Info(proc) % n+1
        Info(proc) % Particles(Info(proc) % n) = i
        Info(proc) % Gindex(Info(proc) % n) = PI % GlobalDOFs(node)
      END DO
    END DO
    
    n = 0
    DO i=1,NoPartitions
      n = n + Info(i) % n
    END DO
    n = 2*(n+2*(2*n*dim+n) + MPI_BSEND_OVERHEAD*2*NoPartitions)
    CALL CheckBuffer(n)
    
    ! Send particles:
    ! ---------------
    DO j=1,NoPartitions
      n = Info(j) % n
      IF ( n<=0 ) CYCLE
      
      CALL MPI_BSEND( Info(j) % Gindex, n, MPI_INTEGER, &
          Neigh(j), 2001, ELMER_COMM_WORLD, ierr )
      
      ALLOCATE(Buf(2*n*dim+n))
      m = 0
      DO k=1,dim
        DO l=1,n
          m = m + 1
          Buf(m) = Particles % Coordinate(Info(j) % Particles(l),k)
        END DO
      END DO
!      DO l=1,n
!        m = m + 1
!        Buf(m) = Particles % Dt(Info(j) % Particles(l))
!      END DO
      IF ( ASSOCIATED(Particles % Velocity) ) THEN
        DO k=1,dim
          DO l=1,n
            m = m + 1
            Buf(m) = Particles % Velocity(Info(j) % Particles(l),k)
          END DO
        END DO
      END IF
      CALL MPI_BSEND( Buf, m, MPI_DOUBLE_PRECISION, &
          Neigh(j), 2002, ELMER_COMM_WORLD, ierr )
      DEALLOCATE(Buf)
    END DO
    
    CALL MPI_WaitAll( NoPartitions, Requests, MPI_STATUSES_IGNORE, ierr )
    n = SUM(Recv_Parts)
    IF ( Particles % NumberOfParticles+n > Particles % MaxNumberOfParticles ) THEN
      CALL IncreaseParticles( Particles, Particles % NumberOfParticles+2*n - &
          Particles % MaxNumberOfParticles )
    END IF
    
    
    ! Recv particles:
    ! ---------------
    DO i=1,NoPartitions
      n = Recv_Parts(i)
      IF ( n<=0 ) CYCLE
      
      proc = Neigh(i)
      
      ALLOCATE(Indexes(n))
      CALL MPI_RECV( Indexes, n, MPI_INTEGER, proc, &
          2001, ELMER_COMM_WORLD, status, ierr )
      
      n_part=Particles % NumberOfParticles
      DO j=1,n
        n_part = n_part+1
        Particles % Status(n_part) = PARTICLE_GHOST
        node = SearchNode(PI,Indexes(j))
        IF ( node<=0 ) STOP 'a'
        Particles % ClosestNode(n_part) = node
        Particles % NoClosestParticle(node) = &
            Particles % NoClosestParticle(node) + 1
      END DO
      DEALLOCATE(Indexes)
      
      m = n+n*dim
      IF ( ASSOCIATED(Particles % Velocity) ) m=m+n*dim
      
      ALLOCATE(Buf(m))
      CALL MPI_RECV( Buf, m, MPI_DOUBLE_PRECISION, proc, &
          2002, ELMER_COMM_WORLD, status, ierr )
      
      n_part=Particles % NumberOfParticles
      m = 0
      DO k=1,dim
        DO l=1,n
          m = m + 1
          Particles % Coordinate(n_part+l,k)=Buf(m)
        END DO
      END DO
            
      IF ( ASSOCIATED(Particles % Velocity) ) THEN
        DO k=1,dim
          DO l=1,n
            m = m + 1
            Particles % Velocity(n_part+l,k)=Buf(m)
          END DO
        END DO
      END IF
      DEALLOCATE(Buf)
      Particles % NumberOfParticles = Particles % NumberOfParticles+n
    END DO
    
    DEALLOCATE(Perm)
    DO i=1,NoPartitions
      IF ( Info(i) % n==0 ) CYCLE
      DEALLOCATE( Info(i) % Gindex, Info(i) % Particles )
    END DO
    DEALLOCATE(Info, Recv_Parts, Neigh, Requests)
    
  END SUBROUTINE CreateGhostParticles
  !------------------------------------------------------------


  !---------------------------------------------------------------
  !> This assumes that the real particles are followed by ghost particles. 
  !> which are destroyed before leaving this configuration.
  !---------------------------------------------------------------
  SUBROUTINE DestroyGhostParticles( Particles ) 

    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No, NumberOfParticles, FirstGhost

    NumberOfParticles = Particles % NumberOfParticles
    FirstGhost = Particles % FirstGhost

    IF ( FirstGhost <= NumberOfParticles ) THEN
      DO No = FirstGhost, NumberOfParticles
        Particles % Status(No) = PARTICLE_LOST
      END DO
      Particles % NumberOfParticles = FirstGhost - 1
    END IF

  END SUBROUTINE DestroyGhostParticles



  !------------------------------------------------------------
  !> For the first call of given node do the list, thereafter 
  !> Return the index until the list is finished.
  !------------------------------------------------------------
  FUNCTION GetNextNeighbour( Particles, No ) RESULT ( No2 )
    IMPLICIT NONE
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER :: No, No2
    
    INTEGER :: PrevNo = 0
    INTEGER, POINTER :: NodeIndexes(:), NeighbourList(:) => NULL(), TmpList(:) => NULL()
    INTEGER :: i,j,k,n,ListSize,NoNeighbours,ElementIndex,Cnt
    LOGICAL :: Visited = .FALSE.
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    
    SAVE Visited,PrevNo,NeighbourList,ListSize,NoNeighbours,Cnt
    
    IF( PrevNo /= No ) THEN
      PrevNo = No
      IF( .NOT. Visited ) THEN
        Visited = .TRUE.
        Mesh => GetMesh()
        n = Mesh % MaxElementNodes 
        ListSize = n * Particles % MaxClosestParticles + 10
        ALLOCATE( NeighbourList( ListSize ) )
        NeighbourList = 0 
        Mesh => GetMesh()
      END IF
      
      Mesh => GetMesh()
      ElementIndex = Particles % ElementIndex(No)
      Element => Mesh % Elements( ElementIndex )
      n = GetElementNOFNodes(Element)
      NodeIndexes => Element % NodeIndexes
      
      NoNeighbours = 0
      DO i=1,n
        j = NodeIndexes(i)
        
        DO k=Particles % CumClosestParticle(j),Particles % CumClosestParticle(j+1)-1 
          No2 = Particles % ClosestParticle(k)
          
          ! No self coupling in this list	
          IF( No2 == No ) CYCLE
          
          ! Set symmetric forces Fij=-Fij so no need to go through twice
          IF ( No2 < No ) CYCLE
          
          NoNeighbours = NoNeighbours + 1
          
          IF( NoNeighbours > ListSize ) THEN
            ALLOCATE( TmpList( ListSize + 20 ) )
            TmpList(1:ListSize) = NeighbourList
            DEALLOCATE( NeighbourList ) 
            NeighbourList => TmpList
            ListSize = ListSize + 20
            NULLIFY( TmpList ) 
            CALL Info('GetNextNeighbour','Allocating more space: '//I2S(ListSize))
          END IF
          
          NeighbourList(NoNeighbours) = No2
        END DO
      END DO
      Cnt = 0
    END IF
    
    Cnt = Cnt + 1
    IF( Cnt > NoNeighbours ) THEN
      No2 = 0
    ELSE
      No2 = NeighbourList( Cnt ) 
    END IF
    
  END FUNCTION GetNextNeighbour
  
  
!------------------------------------------------------------
!> Computes interaction between particles given an interaction 
!> kernel that that is the pointer to the function that computes
!> the effect of the interaction in terms of forces or new coordinate
!> values (for collision models). 
!------------------------------------------------------------ 
  SUBROUTINE ParticleParticleInteraction( Particles, dtime, Collision, InteractionKernel ) 

    IMPLICIT NONE

    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dP) :: dtime
    LOGICAL :: Collision

    INTERFACE 
      SUBROUTINE InteractionKernel( t, r, r2, v, v2, f, f2, Hit ) 
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
        REAL(KIND=dp) :: t,r(3),r2(3),v(3),v2(3),f(3),f2(3)
        LOGICAL :: Hit
      END SUBROUTINE InteractionKernel
    END INTERFACE
    

    INTEGER :: No, No2
    REAL(KIND=dp) :: Coord(3), Velo(3), Coord2(3), Velo2(3), Force(3), Force2(3)
    LOGICAL :: Interact
    
    Coord = 0.0_dp
    Velo = 0.0_dp
    Force = 0.0_dp    
    Coord2 = 0.0_dp
    Velo2 = 0.0_dp
    Force2 = 0.0_dp
    

    DO No=1,Particles % NumberOfParticles
      IF ( Particles % Status(no) == PARTICLE_GHOST ) EXIT
      IF ( Particles % Status(no) == PARTICLE_LOST  ) CYCLE
      
      Coord = GetParticleCoord( Particles, No )
      Velo  = GetParticleVelo( Particles, No )

      DO WHILE(.TRUE.)
        No2 = GetNextNeighbour( Particles, No )           
        IF( No2 == 0 ) EXIT
        Coord2 = GetParticleCoord( Particles, No2 )
        Velo2  = GetParticleVelo( Particles, No2 )
        
        CALL InteractionKernel(dtime,Coord,Coord2,Velo,Velo2,&
            Force,Force2, Interact) 
        IF(.NOT. Interact ) CYCLE
          
        IF( Collision ) THEN
          CALL SetParticleCoord( Particles, No, Coord )
          CALL SetParticleCoord( Particles, No2, Coord2 )
          CALL AddParticleForce( Particles, No, Force )
          CALL AddParticleForce( Particles, No2, Force2 )
        ELSE
          CALL AddParticleForce( Particles, No, Force )
          CALL AddParticleForce( Particles, No2, Force2 )
        END IF

      END DO
    END DO
  END SUBROUTINE ParticleParticleInteraction


  !---------------------------------------------------------
  !> Initialize the time for next time integration step
  !---------------------------------------------------------
  SUBROUTINE ParticleInitializeTime( Particles, No )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: No
    
    IF( PRESENT( No ) ) THEN
      Particles % Force( No, : ) = 0.0_dp
    ELSE
      Particles % Force = 0.0_dp
    END IF
    
  END SUBROUTINE ParticleInitializeTime

 

  !---------------------------------------------------------
  !> Advance the particles with a time step. The timestep may
  !> also be an intermediate Runge-Kutta step.
  !---------------------------------------------------------
  SUBROUTINE ParticleAdvanceTimestep( Particles, RKstepInput )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: RKStepInput

    REAL(KIND=dp) :: dtime
    TYPE(Variable_t), POINTER :: Var, TimeVar, DistVar, DtVar
    LOGICAL :: GotVar, GotTimeVar, GotDistVar, MovingMesh
    REAL(KIND=dp) :: ds, dCoord(3),Coord(3),Velo(3),Speed0,Speed
    INTEGER :: dim, Status, TimeOrder, No, NoMoving
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoParticles, RKStep
    LOGICAL :: Found, Visited = .FALSE.,RK2,HaveSpeed0

    REAL(KIND=dp) :: mass, drag
    REAL(KIND=dp), POINTER :: massv(:), dragv(:)
    LOGICAL :: GotMass, GotDrag
    INTEGER :: CurrGroup, PrevGroup, NoGroups
    CHARACTER(*), PARAMETER :: Caller = 'ParticleAdvanceTimestep'

    SAVE TimeOrder, dim, Mass, Drag, Visited, dCoord, Coord, GotTimeVar, &
	GotDistVar, TimeVar, DtVar, DistVar, MovingMesh,Speed0,HaveSpeed0, Params

    
    IF(.NOT. Visited ) THEN
      Params => ListGetSolverParams()
      TimeOrder = Particles % TimeOrder
      dim = Particles % dim
      
      dCoord = 0.0_dp
      Coord = 0.0_dp
      Visited = .TRUE.
      MovingMesh = .FALSE.

      TimeVar => ParticleVariableGet( Particles,'particle time')
      GotTimeVar = ASSOCIATED( TimeVar )

      IF( GotTimeVar ) THEN
        IF( .NOT. Particles % DtConstant ) THEN
          DtVar => ParticleVariableGet( Particles,'particle dt')
          IF(.NOT. ASSOCIATED( DtVar ) ) THEN
            CALL Fatal(Caller,'Variable timestep, > particle dt < should exist!')
          END IF
        END IF        
      END IF
 
      DistVar => ParticleVariableGet( Particles,'particle distance')
      GotDistVar = ASSOCIATED( DistVar )

      Speed0 = GetCReal( Params,'Particle Speed Constant',HaveSpeed0 )      
    END IF

    NoParticles = Particles % NumberOfParticles
    NoGroups = Particles % NumberOfGroups     
    NoMoving = 0
    RK2 = Particles % RK2
    RKStep = 0
    IF(PRESENT(RKStepInput)) RKStep=RKStepInput

    IF( RK2 .AND. .NOT. ASSOCIATED( Particles % PrevVelocity ) ) THEN
      ALLOCATE( Particles % PrevVelocity( &
          SIZE( Particles % Velocity,1 ),SIZE( Particles % Velocity,2) ) )
      Particles % PrevVelocity = Particles % Velocity
    END IF

    IF( Particles % DtConstant ) THEN
      dtime = Particles % DtSign * Particles % dTime
    END IF

    GotMass = .FALSE.
    GotDrag = .FALSE.
    IF( TimeOrder == 2 ) THEN      
      IF( NoGroups > 1 ) THEN
        massv => ListGetConstRealArray1( Params,'Particle Mass',GotMass)
        Mass = 0.0_dp
      ELSE
        Mass = ListGetConstReal( Params,'Particle Mass',GotMass)
      END IF
      IF(.NOT. GotMass) CALL Fatal(Caller,'> Particle Mass < should be given!')
    ELSE IF( TimeOrder == 1 ) THEN
      IF( NoGroups > 1 ) THEN
        dragv => ListGetConstRealArray1( Params,'Particle Drag Coefficient',GotDrag)
        Drag = 0.0_dp
      ELSE
        Drag = ListGetConstReal( Params,'Particle Drag Coefficient',GotDrag)
      END IF
      IF(.NOT. GotDrag) CALL Fatal(Caller,'> Particle Drag Coefficient < should be given!')
    END IF
    PrevGroup = -1 

    
    
    ! Now move the particles
    !---------------------------
    DO No=1, NoParticles

      Status = Particles % Status(No)

      IF ( Status >= PARTICLE_LOST ) CYCLE
      IF ( Status <= PARTICLE_INITIATED ) CYCLE
      IF ( Status == PARTICLE_WALLBOUNDARY ) CYCLE

      ! Cumulate the time
      !-----------------------------------------------------------------       
      IF( .NOT. Particles % DtConstant ) THEN
        dtime = Particles % DtSign * DtVar % Values(No)
        TimeVar % Values(No) = TimeVar % Values(No) + DtVar % Values(No)
        IF( ABS( dtime ) < TINY( dtime ) ) CYCLE
      ELSE IF( GotTimeVar ) THEN
        TimeVar % Values(No) = TimeVar % Values(No) + Particles % dTime         
      END IF

      IF( NoGroups > 1 ) THEN
        CurrGroup = GetParticleGroup(Particles,No)
        IF( CurrGroup /= PrevGroup ) THEN
          IF(GotMass) THEN
            mass = massv(MIN(SIZE(massv),CurrGroup)) 
          END IF
          IF(GotDrag) drag = dragv(MIN(SIZE(dragv),CurrGroup))
        END IF
        PrevGroup = CurrGroup
      END IF
      
      IF ( Status == PARTICLE_FIXEDCOORD ) THEN
        Particles % Velocity(No,:) = 0.0_dp
	CYCLE
      ELSE IF( Status == PARTICLE_FIXEDVELO ) THEN
        CONTINUE
      ELSE IF( TimeOrder == 2 ) THEN
        Particles % Velocity(No,:) = Particles % Velocity(No,:) + &
            dtime * Particles % Force(No,:) / Mass
      ELSE IF( TimeOrder == 1 ) THEN
        Particles % Velocity(No,:) = Particles % Force(No,:) / Drag      
      ELSE IF( TimeOrder == 0 ) THEN
        ! Velocity stays fixed
        CONTINUE
      ELSE
        CALL Fatal(Caller,'Unknown time order')
      END IF
       
         
      IF( RK2 .AND. RKStep == 2 ) THEN
         Velo(1:dim) = &
           ( 2 * Particles % Velocity(No,:) - Particles % PrevVelocity(No,:) )
       ELSE
        Velo(1:dim) = Particles % Velocity(No,:) 	
      END IF

      IF( HaveSpeed0 ) THEN
        Speed = SQRT( SUM( Velo(1:dim)**2 ) )
        IF( Speed > TINY( Speed ) ) THEN
          dCoord(1:dim) = dtime * Speed0 * Velo(1:dim) / Speed
        ELSE
          dCoord(1:dim) = 0.0_dp
        END IF
      ELSE
        dCoord(1:dim) = dtime * Velo(1:dim) 
      END IF

      Particles % PrevCoordinate(No,:) = Particles % Coordinate(No,:)
      IF( ASSOCIATED( Particles % PrevVelocity ) ) THEN
        Particles % PrevVelocity(No,:) = Particles % Velocity(No,:)
      END IF
      Particles % Force(No,:)= 0.0_dp
 
      NoMoving = NoMoving + 1

      Particles % Status(No) = PARTICLE_READY
      Particles % Coordinate(No,:) = Particles % Coordinate(No,:) + dCoord(1:dim)

      IF( GotDistVar ) THEN
        ds = SQRT( SUM( dCoord(1:dim)**2 ) )
        DistVar % Values(No) = DistVar % Values(No) + ds    
      END IF
    END DO

    IF( Particles % DtConstant ) THEN
      Particles % Time = Particles % Time + Particles % dTime
    END IF
    
    Particles % NumberOfMovingParticles = NoMoving

 END SUBROUTINE ParticleAdvanceTimestep


  !---------------------------------------------------------
  !> Advance some tracer quantities related to the particles.
  !---------------------------------------------------------
  SUBROUTINE ParticlePathIntegral( Particles, RKstepInput )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: RKstepInput

    TYPE(Variable_t), POINTER :: TimeIntegVar, DistIntegVar, DtVar, &
        PartTimeVar, MeshDtVar 
    LOGICAL :: GotVar, RK2
    REAL(KIND=dp) :: ds,dtime,Coord(3),PrevCoord(3),LocalCoord(3),Velo(3),u,v,w,&
        Source,PrevSource,detJ,RKCoeff,GradSource(3),MeshDt, dtRat, DummyVals(1)
    INTEGER :: dim, Status, RKStep
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoParticles, No, n, NoVar, i, j, bf_id
    LOGICAL :: Found, Stat, Visited = .FALSE., UseDummy
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    INTEGER, POINTER :: Indexes(:)
    TYPE(ValueList_t), POINTER :: BodyForce
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: TimeInteg = .FALSE., DistInteg = .FALSE., UseGradSource = .FALSE., TimeDepFields = .FALSE.
    TYPE(ValueHandle_t), SAVE :: TimeSource_h, DistSource_h
    CHARACTER(:), ALLOCATABLE :: VariableName
    CHARACTER(*), PARAMETER :: Caller = 'ParticlePathIntegral'


    SAVE TimeInteg, DistInteg, dim, Visited, Mesh, DtVar, Basis, Nodes, Params, &
        TimeIntegVar, DistIntegVar, UseGradSource, dBasisdx, TimeDepFields, &
        PartTimeVar, MeshDtVar, UseDummy

    CALL Info(Caller,'Integrating variables over the path',Level=12)


    ! If Runge-Kutta is used take the mid-point rule.
    RKSTep = 0
    IF( PRESENT( RKStepInput ) ) RKStep = RKStepInput
    
    IF( RKStep > 1 ) RETURN
    
    IF(.NOT. Visited ) THEN
      Visited = .TRUE.

      TimeIntegVar => ParticleVariableGet( Particles,'particle time integral')
      TimeInteg = ASSOCIATED( TimeIntegVar )
      IF( TimeInteg ) THEN        
        IF( .NOT. ListCheckPresentAnyBodyForce( CurrentModel,&
            'Particle Time Integral Source') ) THEN
          CALL Fatal(Caller,'Path integral requires body force: "Particle Time Integral Source"')
        END IF
      END IF
      
      DistIntegVar => ParticleVariableGet( Particles,'particle distance integral')
      DistInteg = ASSOCIATED( DistIntegVar )
      IF( DistInteg ) THEN
        IF( .NOT. ListCheckPresentAnyBodyForce( CurrentModel,&
            'Particle Distance Integral Source') ) THEN
          CALL Fatal(Caller,'Path integral requires body force: "Particle Distance Integral Source"')
        END IF
      END IF

      IF( .NOT. (TimeInteg .OR. DistInteg ) ) RETURN
      
      Params => ListGetSolverParams()
      Mesh => CurrentModel % Solver % Mesh
      dim = Particles % dim

      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n), dBasisdx(n,3) )
      Basis = 0.0_dp
      Nodes % x = 0.0_dp
      Nodes % y = 0.0_dp
      Nodes % z = 0.0_dp
      dBasisdx = 0.0_dp
      
      UseGradSource = ListGetLogical( Params,'Source Gradient Correction',Found)
      ! If the correction is not given follow the logic of velocity estimation
      IF(UseGradSource .AND. Particles % RK2 ) THEN
        CALL Fatal(Caller,'Quadratic source correction incompatibe with Runge-Kutta')
      END IF

      TimeDepFields = ListGetLogical( Params,'Source Time Correction',Found ) 
      IF(TimeDepFields ) THEN
        IF( Particles % RK2 ) THEN
          CALL Fatal(Caller,'Time correction is incompatibe with Runge-Kutta')
        END IF
        PartTimeVar => ParticleVariableGet( Particles, 'particle time' )
        MeshDtVar => VariableGet( Mesh % Variables,'timestep size')
      END IF
        
      IF( .NOT. Particles % DtConstant ) THEN
        DtVar => ParticleVariableGet( Particles,'particle dt')
        IF(.NOT. ASSOCIATED( DtVar ) ) THEN
          CALL Fatal(Caller,'Variable timestep, > particle dt < should exist!')
        END IF
      END IF

      UseDummy = ListGetLogical( Params,'Particle Integral Dummy Argument',Found )

      IF( UseDummy .AND. UseGradSource ) THEN
        CALL Fatal(Caller,'Gradient correction and dummy arguments are incompatible for now!')
      END IF
    END IF
      
    ! Nothing to integrate over
    IF( .NOT. (TimeInteg .OR. DistInteg ) ) RETURN

    j = 0
    IF( UseDummy ) THEN
      CALL Info(Caller,'Expecting one dummy argument for integral sources!',Level=30)
      j = 1
    END IF
      
    IF(TimeInteg) THEN
      CALL ListInitElementKeyword( TimeSource_h,'Body Force','Particle Time Integral Source',DummyCount=j)    
    ELSE
      CALL ListInitElementKeyword( DistSource_h,'Body Force','Particle Distance Integral Source',DummyCount=j)
    END IF
      
    NoParticles = Particles % NumberOfParticles
    RK2 = Particles % RK2
    IF( RK2 ) THEN
      RKCoeff = 2.0
    ELSE
      RKCoeff = 1.0
    END IF

    IF( dim < 3 ) w = 0.0_dp

    IF( Particles % DtConstant ) THEN
      dtime = RKCoeff * Particles % dTime
    END IF
    MeshDt = 0.0_dp
    IF( TimeDepFields ) THEN
      MeshDt = MeshDtVar % Values(1)
    END IF

    j = 0    
    DO No=1, NoParticles
      Status = Particles % Status(No)
      
      IF ( Status >= PARTICLE_LOST ) CYCLE
      IF ( Status <= PARTICLE_INITIATED ) CYCLE
      IF ( Status == PARTICLE_WALLBOUNDARY ) CYCLE
      IF ( Status == PARTICLE_FIXEDCOORD ) CYCLE

      ! Local timestep size
      !-----------------------------------------------------------------       
      IF( .NOT. Particles % DtConstant ) THEN
        dtime = RKCoeff * DtVar % Values(No)
      END IF

      ! This is the ratio of time vs. the mesh timestep. Starts from 0 and goes to 1. 
      IF( TimeDepFields ) THEN
        dtrat = PartTimeVar % Values(No) / MeshDt
      END IF
      
      Coord = 0._dp
      Coord(1:dim) = Particles % Coordinate(No,:) 
      Velo  = 0._dp
      Velo(1:dim) = Particles % Velocity(No,:) 	
      
      Element => Mesh % Elements( Particles % ElementIndex(No) )            
      CurrentModel % CurrentElement => Element

      n = Element % TYPE % NumberOfNodes
      Indexes => Element % NodeIndexes

      Nodes % x(1:n) = Mesh % Nodes % x( Indexes ) 
      Nodes % y(1:n) = Mesh % Nodes % y( Indexes ) 
      Nodes % z(1:n) = Mesh % Nodes % z( Indexes ) 

      bf_id = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values, &
          'Body Force', Found )
      IF( .NOT. Found ) CYCLE   
      BodyForce => CurrentModel % BodyForces(bf_id) % Values

      IF( ASSOCIATED( Particles % uvw ) ) THEN
        u = Particles % uvw(No,1)
        v = Particles % uvw(No,2)
        IF( dim == 3 ) w = Particles % uvw(No,3)
      ELSE
        IF(.NOT. PointInElement( Element, Nodes, &
            Coord, LocalCoord ) ) CYCLE
        u = LocalCoord(1)
        v = LocalCoord(2)
        w = LocalCoord(3)
      END IF

      ! Make the integral more accurate by including a correction term based on the 
      ! elemental derivate of the source term. 
      IF( UseGradSource ) THEN
        stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis,dBasisdx)
      ELSE
        stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis)
      END IF

      PrevCoord(1:dim) = Particles % PrevCoordinate(No,:) 

      ! Path integral over time
      IF( TimeInteg ) THEN        
        IF( UseDummy ) THEN
          DummyVals(1) = TimeIntegVar % Values(No)
          Source = ListGetElementReal( TimeSource_h, Basis, Element, Found, DummyVals = DummyVals )         
        ELSE
          Source = ListGetElementReal( TimeSource_h, Basis, Element, Found )
        END IF
        IF( Found ) THEN
          IF ( UseGradSource ) THEN
            GradSource = ListGetElementRealGrad( TimeSource_h,dBasisdx,Element)                    
            Source = Source + 0.5*SUM( GradSource(1:dim) * (PrevCoord(1:dim) - Coord(1:dim)) )
          END IF

          IF( TimeDepFields ) THEN
            IF( UseDummy ) THEN
              PrevSource = ListGetElementReal( TimeSource_h, Basis, Element, Found,tstep=-1, DummyVals = DummyVals  )
            ELSE
              PrevSource = ListGetElementReal( TimeSource_h, Basis, Element, Found,tstep=-1 )
            END IF
            IF ( UseGradSource ) THEN
              GradSource = ListGetElementRealGrad( TimeSource_h,dBasisdx,Element,tstep=-1)                    
              PrevSource = PrevSource + 0.5*SUM( GradSource(1:dim) * (PrevCoord(1:dim) - Coord(1:dim)) )
            END IF
            TimeIntegVar % Values(No) = TimeIntegVar % Values(No) + dtime * &
                ( (1-dtrat)*Source + dtrat * PrevSource ) 
          ELSE          
            TimeIntegVar % Values(No) = TimeIntegVar % Values(No) + dtime * Source
          END IF
        END IF        
      END IF
      
      ! Path integral over distance
      IF( DistInteg ) THEN
        IF( RK2 ) THEN
          ! for R-K the velocity has been updated to the midpoint and this is used
          ! to determine the differential path. 
          ds = dtime * SQRT( SUM( Velo(1:dim)**2) )
        ELSE
          ! If we have used quadtatic velocity correction at the previous point
          ! the current (or previous) velocity alone does not give the correct 
          ! ds, but this does.
          ds = SQRT(SUM((PrevCoord(1:dim) - Coord(1:dim))**2))
        END IF
        
        IF( UseDummy ) THEN
          DummyVals(1) = DistIntegVar % Values(No)
          Source = ListGetElementReal( DistSource_h, Basis, Element, Found, DummyVals = DummyVals )         
        ELSE
          Source = ListGetElementReal( DistSource_h, Basis, Element, Found )
        END IF

        IF( Found ) THEN
          IF ( UseGradSource ) THEN
            GradSource = ListGetElementRealGrad( DistSource_h,dBasisdx,Element)                  
            Source = Source + 0.5*SUM( GradSource(1:dim) * (PrevCoord(1:dim) - Coord(1:dim)) )
          END IF          
          IF( TimeDepFields ) THEN
            IF( UseDummy ) THEN
              PrevSource = ListGetElementReal( DistSource_h, Basis, Element, Found,tstep=-1, DummyVals = DummyVals  )
            ELSE
              PrevSource = ListGetElementReal( DistSource_h, Basis, Element, Found,tstep=-1  )
            END IF
            IF ( UseGradSource ) THEN
              GradSource = ListGetElementRealGrad( DistSource_h,dBasisdx,Element,tstep=-1)                    
              PrevSource = PrevSource + 0.5*SUM( GradSource(1:dim) * (PrevCoord(1:dim) - Coord(1:dim)) )
            END IF
            DistIntegVar % Values(No) = DistIntegVar % Values(No) + ds * &
                ( (1-dtrat)*Source + dtrat * PrevSource ) 
          ELSE          
            DistIntegVar % Values(No) = DistIntegVar % Values(No) + ds * Source
          END IF
        END IF
      END IF

    END DO

  END SUBROUTINE ParticlePathIntegral



  !---------------------------------------------------------------    
  !> Checks the boundaries for rectangular and hexahedral shapes and 
  !> enforces periodic BCs. Currently the only supported way
  !> for setting periodic BCs.
  !---------------------------------------------------------------    
  SUBROUTINE ParticleBoxPeriodic( Particles )
    
    TYPE(Particle_t), POINTER :: Particles
    
    TYPE(Solver_t), POINTER :: Solver
    REAL(KIND=dp) :: Coord, Rad
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    REAL(KIND=dP) :: MinCoord(3), MaxCoord(3), EpsCoord
    INTEGER :: i,j,k,dim, ierr, PeriodicDir(3),NoPeriodic
    LOGICAL :: Mapped,Reflect,Found,SaveCount,Visited = .FALSE.
    INTEGER :: Operations, No, NoParticles, Status, NoCount(6), NoStep
    INTEGER, POINTER :: TmpInteger(:)
    CHARACTER(:), ALLOCATABLE :: Filename
    
    SAVE Visited, Reflect, PeriodicDir, NoPeriodic, MinCoord, MaxCoord, dim, &
        SaveCount, NoCount, Filename, NoStep
    
    IF( .NOT. Visited ) THEN
      Visited = .TRUE.
      Mesh => GetMesh()
      Params => ListGetSolverParams()
      dim = Mesh % Meshdim
      
      NoPeriodic = 0
      PeriodicDir = 0
      
      TmpInteger => ListGetIntegerArray( &
          Params,'Box Periodic Directions',Found )     
      IF( Found ) THEN
        NoPeriodic = SIZE( TmpInteger )
        DO i=1,NoPeriodic
          PeriodicDir(i) = TmpInteger(i)
        END DO
      ELSE IF( ListGetLogical( Params,'Box Particle Periodic',Found)) THEN
        NoPeriodic = dim
        DO i=1,dim
          PeriodicDir(i) = i
        END DO
      END IF
      
      MinCoord = Particles % GlobalMinCoord 
      MaxCoord = Particles % GlobalMaxCoord 
      EpsCoord = EPSILON( EpsCoord ) * MAXVAL( MaxCoord - MinCoord )
      MinCoord = MinCoord + EpsCoord 
      MaxCoord = MaxCoord - EpsCoord

      Filename = ListGetString( Params,'Box Periodic Filename', SaveCount )
      NoCount = 0
      NoStep = 0
    END IF
    
    IF( NoPeriodic == 0 ) RETURN
    
    NoParticles = Particles % NumberOfParticles
    
    
    DO No = 1, NoParticles
      Status = Particles % Status(No) 
      IF( Status >= PARTICLE_LOST ) CYCLE
      IF( Status < PARTICLE_INITIATED ) CYCLE
      
      ! Boundary conditions for periodic BCs
      !------------------------------------------
      DO i=1,NoPeriodic
        Mapped = .FALSE.
        DO j=1,NoPeriodic
          k = PeriodicDir(j)
          coord = Particles % Coordinate(No,k)
          IF( coord < MinCoord(k) ) THEN
            IF( SaveCount ) NoCount(2*k-1) = NoCount(2*k-1) + 1
            coord = MaxCoord(k) - MinCoord(k) + coord
            Particles % Coordinate(No,k) = coord
            Mapped = .TRUE.
          ELSE IF ( coord > MaxCoord(k) ) THEN
            IF( SaveCount ) NoCount(2*k) = NoCount(2*k) + 1 
            Coord = MinCoord(k) - MaxCoord(k) + Coord
            Particles % Coordinate(No,k) = coord           
            Mapped = .TRUE.
          END IF
        END DO
        
        IF(.NOT. Mapped ) EXIT
      END DO
    END DO
    
    IF( SaveCount ) THEN      
      IF( NoStep == 0 ) THEN
        OPEN (10, FILE=FileName )
      ELSE
        OPEN (10, FILE=FileName, POSITION='append')
      END IF
      NoStep = NoStep + 1
      WRITE( 10, * ) NoStep, NoCount(1:2*NoPeriodic)
      CLOSE( 10 )
    END IF

  END SUBROUTINE ParticleBoxPeriodic
  
  
  !---------------------------------------------------------------    
  !> Checks the boundaries for rectangular and hexahedral shapes and 
  !> enforces elastic reflection. This is alternative and 
  !> computationally more economic way and is ideal for testing
  !> purposes, at least.
  !---------------------------------------------------------------    
  SUBROUTINE ParticleBoxContact(Particles)
    
    TYPE(Particle_t), POINTER :: Particles
    
    TYPE(Solver_t), POINTER :: Solver
    REAL(KIND=dp) :: Coord, Velo, Rad, Spring
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    REAL(KIND=dP) :: MinCoord(3), MaxCoord(3), eta
    INTEGER :: i,j,k,dim, ierr,ContactDir(3)
    LOGICAL :: Mapped,Found,CollisionBC,ContactBC,Visited = .FALSE.
    INTEGER :: No, NoParticles, Status, NoContact
    INTEGER, POINTER :: TmpInteger(:)
    
    SAVE Visited, NoContact, ContactDir, MinCoord, MaxCoord, dim, &
        CollisionBC, ContactBC, Spring
    
    IF( .NOT. Visited ) THEN
      Visited = .TRUE.
      Mesh => GetMesh()
      Params => ListGetSolverParams()
      dim = Mesh % Meshdim
      
      NoContact = 0
      ContactDir = 0
      
      ContactBC = ListGetLogical(Params,'Box Particle Contact',Found)
      CollisionBC = ListGetLogical( Params,'Box Particle Collision',Found)
      
      IF( ContactBC .OR. CollisionBC ) THEN
        TmpInteger => ListGetIntegerArray( &
            Params,'Box Contact Directions',Found )     
        IF( Found ) THEN
          NoContact = SIZE( TmpInteger )
          DO i=1,NoContact
            ContactDir(i) = TmpInteger(i)
          END DO
        ELSE
          DO i=1,dim
            ContactDir(i) = i
          END DO
          NoContact = dim
        END IF
      ELSE
        NoContact = 0
      END IF
      IF( NoContact == 0 ) RETURN
      
      MinCoord = Particles % GlobalMinCoord
      MaxCoord = Particles % GlobalMaxCoord
      
      ! Particles of finite size collide before their center 
      ! hits the wall.
      Rad = GetCReal( Params,'Wall Particle Radius',Found)    
      IF( Found ) THEN
        MaxCoord = MaxCoord - Rad
        MinCoord = MinCoord + Rad
      END IF
      
      IF( ContactBC ) THEN
        Spring = GetCReal(Params,'Wall Particle Spring',Found)
        IF(.NOT. Found) CALL Fatal('ParticleBoxContact',&
            '> Wall Particle Spring < needed!')
      END IF
      
    END IF
    
    IF( NoContact == 0) RETURN
    
    NoParticles = Particles % NumberOfParticles
    
    DO No = 1, NoParticles
      Status = Particles % Status(No) 
      IF( Status >= PARTICLE_LOST ) CYCLE
      IF( Status < PARTICLE_INITIATED ) CYCLE
      
      ! Boundary conditions for reflection. 
      ! Multiple reflections may be carried out.
      !------------------------------------------
      DO i=1,NoContact
        
        IF( CollisionBC ) THEN        
          Mapped = .FALSE.
          DO j=1,NoContact
            k = ContactDir(j)
            Coord = Particles % Coordinate(No,k)
            
            IF( Coord < MinCoord(k) ) THEN
              Coord = 2 * MinCoord(k) - Coord
              Particles % Coordinate(No,k) = Coord
              Particles % Velocity(No,k) = -Particles % Velocity(No,k)
              Mapped = .TRUE.
            ELSE IF ( Coord > MaxCoord(k) ) THEN
              Coord = 2 * MaxCoord(k) - Coord
              Particles % Coordinate(No,k) = Coord          
              Particles % Velocity(No,k) = -Particles % Velocity(No,k)
              Mapped = .TRUE.
            END IF
          END DO
          IF(.NOT. Mapped ) EXIT
        ELSE        
          k = ContactDir(i)
          Coord = Particles % Coordinate(No,k)
          
          IF( MinCoord(k) - Coord > 0.0_dp ) THEN
            eta = MinCoord(k) - Coord          
            Particles % Force(No,k) = Particles % Force(No,k) + eta * Spring
          ELSE IF( Coord - MaxCoord(k) > 0.0_dp ) THEN
            eta = Coord - MaxCoord(k) 
            Particles % Force(No,k) = Particles % Force(No,k) - eta * Spring
          END IF
        END IF
      END DO
      
    END DO
    
  END SUBROUTINE ParticleBoxContact



!--------------------------------------------------------------------------
!> Set a timestep for the particles.
!> Depending on the definitions the timestep may be the same for all 
!> particles, or may be defined independently for each particle.
!-------------------------------------------------------------------------
  FUNCTION GetParticleTimeStep(Particles, InitInterval, tinit ) RESULT ( dtout )
    
    TYPE(Particle_t), POINTER :: Particles
    LOGICAL :: InitInterval
    REAL(KIND=dp), OPTIONAL :: tinit
    REAL(KIND=dp) :: dtout

    INTEGER :: No, Status    
    REAL(KIND=dp) :: dt,dt0,tfin,tprev,dsgoal,hgoal,dtmax,dtmin,dtup,dtlow, &
        CharSpeed, CharTime, dtave, dtmax2, dtmin2
    LOGICAL :: GotIt,TfinIs,NStepIs,DsGoalIs,HgoalIs,HgoalIsUniso,DtIs
    INTEGER :: nstep, TimeStep, PrevTimeStep = -1, nset
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: TimeVar, DtVar
    
    SAVE dt0,dsgoal,hgoal,dtmax,dtmin,DtIs,Nstep,&
        tprev,Tfin,TfinIs,DsGoalIs,HgoalIs,HgoalIsUniso,PrevTimeStep, &
	DtVar,TimeVar


    CALL Info('GetParticleTimestep','Setting timesteps for particles!',Level=20)
    
    dtout = 0.0_dp; dtave = 0._dp; nset = 0

    IF( InitInterval ) THEN
      Params => ListGetSolverParams()
      
      ! directly defined timestep
      dt0 = GetCReal(Params,'Timestep Size',DtIs)
      
      ! Constraint by absolute step size taken (in length units)
      dsgoal = GetCReal( Params,'Timestep Distance',DsGoalIs)
      
      ! Constraint by relative step size taken (1 means size of the element)
      hgoal = GetCReal( Params,'Timestep Courant Number',HGoalIs)

      ! Constraint by relative step size taken (each cartesian direction of element)
      HGoalIsUniso = .FALSE.
      IF(.NOT. HGoalIs ) THEN
        hgoal = GetCReal( Params,'Timestep Unisotropic Courant Number',HGoalIsUniso)
      END IF
      
      Nstep = GetInteger( Params,'Max Timestep Intervals',GotIt)
      IF(.NOT. GotIt) Nstep = 1
      
      ! Constraint timestep directly
      dtmax = GetCReal( Params,'Max Timestep Size',GotIt)
      IF(.NOT. GotIt ) dtmax = HUGE( dtmax ) 
      dtmin = GetCReal( Params,'Min Timestep Size',GotIt)
      IF(.NOT. GotIt ) dtmin = 0.0 
      
      TfinIs = .FALSE.
      IF( GetLogical(Params,'Simulation Timestep Sizes',GotIt) ) THEN
        tfin = GetTimeStepsize()
        TfinIs = .TRUE.
      ELSE
        tfin = GetCReal(Params,'Max Cumulative Time',TfinIs)
      END IF
      
      IF( .NOT. Particles % DtConstant ) THEN
        DtVar => ParticleVariableGet( Particles,'particle dt')       
        IF( .NOT. ASSOCIATED( DtVar ) ) THEN
          CALL ParticleVariableCreate( Particles,'particle dt')
          DtVar => ParticleVariableGet( Particles,'particle dt')       
        END IF

        TimeVar => ParticleVariableGet( Particles,'particle time')       
        IF( .NOT. ASSOCIATED( TimeVar ) ) THEN
          CALL Fatal('GetParticleTimestep','Variable > Particle time < does not exist!')
        END IF
      END IF

      tprev = 0.0_dp
    END IF

    
    ! Get upper and lower constraints for timestep size
    ! These generally depend on the velocity field and mesh
    !--------------------------------------------------------------------
    IF( Particles % DtConstant ) THEN
      IF( DtIs ) THEN
        dt = dt0 
      ELSE IF( DsGoalIs ) THEN
        CharSpeed = CharacteristicSpeed( Particles )
        dt = dsgoal / CharSpeed
      ELSE IF( HgoalIs ) THEN
        CharTime = CharacteristicElementTime( Particles )
        dt = Hgoal * CharTime ! ElementH / Speed
      ELSE IF( tfinIs ) THEN
        dt = tfin / Nstep
      ELSE IF( HgoalIsUniso ) THEN
        CALL Fatal('GetParticleTimesStep','Cannot use unisotropic courant number with constant dt!')
      ELSE
        CALL Fatal('GetParticleTimeStep','Cannot determine timestep size!')
      END IF

      ! Constrain the timestep
      !------------------------------------------------------------------
!     dt = MAX( MIN( dt, dtmax ), dtmin )

      ! Do not exceed the total integration time
      !-------------------------------------------
      IF( PRESENT( tinit ) ) tprev = tinit
      IF( TfinIs .AND. dt + tprev > tfin ) THEN
        dt = tfin - tprev
      END IF
      tprev = tprev + dt
      Particles % dtime = dt
      dtout = dt

      WRITE(Message,'(A,ES12.3)') 'Constant particle timestep:',dtout
      CALL Info('GetParticleTimestep', Message,Level=12)           
    ELSE 
      DtVar % Values = 0.0_dp
      dtave = 0.0_dp
      dtmax2 = -HUGE(dtmax2)
      dtmin2 = HUGE(dtmin2)
      
      DO No = 1, Particles % NumberOfParticles

        Status = Particles % Status( No )
        IF ( Status >= PARTICLE_LOST ) CYCLE
        IF ( Status <= PARTICLE_INITIATED ) CYCLE
        IF ( Status == PARTICLE_WALLBOUNDARY ) CYCLE
        IF ( Status == PARTICLE_FIXEDCOORD ) CYCLE

	tprev = TimeVar % Values(No)        

        IF( DtIs ) THEN
          dt = dt0 
        ELSE IF( DsGoalIs ) THEN
          CharSpeed = CharacteristicSpeed( Particles, No )     
          dt = dsgoal / CharSpeed
        ELSE IF( HgoalIs ) THEN
          CharTime = CharacteristicElementTime( Particles, No )     
          dt = Hgoal * CharTime ! ElementH / Speed
        ELSE IF( tfinIs ) THEN
          dt = tfin / Nstep
        ELSE IF( HgoalIsUniso ) THEN
          CharTime = CharacteristicUnisoTime( Particles, No )     
          dt = Hgoal * CharTime ! ElementH / Speed
        ELSE
          CALL Fatal('GetParticlesTimeStep','Cannot determine timestep size!')
        END IF

        ! Constrain the timestep
        !------------------------------------------------------------------
        dt = MAX( MIN( dt, dtmax ), dtmin )

        ! Do not exceed the total integration time
        !-------------------------------------------
        IF( PRESENT( tinit ) ) tprev = tinit

        IF( TfinIs .AND. dt + tprev > tfin ) THEN
          dt = tfin - tprev
        END IF
        DtVar % Values(No) = dt
        nset = nset + 1
        
        dtave = dtave + dt
        dtmin2 = MIN(dtmin2, dt)
        dtmax2 = MAX(dtmax2, dt)          
      END DO

      nset = ParallelReduction(nset) 
      WRITE(Message,'(A,I0)') 'Timestep set for particles: ',nset
      CALL Info('GetParticleTimestep', Message,Level=12)           

      IF( nset == 0 ) THEN
        ! If no particles are set then the indicative forward timestep becomes zero!
        dtout = 0.0_dp
      ELSE        
        dtmax2 = ParallelReduction(dtmax2,2) 
        IF( InfoActive(12) ) THEN
          dtave = ParallelReduction(dtave) / nset 
          dtmin2 = ParallelReduction(dtmin2,1) 
                    
          WRITE(Message,'(A,ES12.3)') 'Average particle timestep:',dtave
          CALL Info('GetParticleTimestep', Message)

          WRITE(Message,'(A,ES12.3)') 'Minimum particle timestep:',dtmin2
          CALL Info('GetParticleTimestep', Message)
        
          WRITE(Message,'(A,ES12.3)') 'Maximum particle timestep:',dtmax2
          CALL Info('GetParticleTimestep', Message)
        END IF
        dtout = dtmax2
      END IF
    END IF
          
    IF( Particles % Rk2 ) THEN
      IF( Particles % DtConstant ) THEN
        Particles % Dtime = 0.5_dp * Particles % Dtime
      ELSE
        DtVar % Values = 0.5_dp * DtVar % Values
      END IF
    END IF
    
  END FUNCTION GetParticleTimeStep



!------------------------------------------------------------------------------
!> Creates a variable related to the particles. The normal variabletype without
!> the permutation vector is used. 
!------------------------------------------------------------------------------
  SUBROUTINE ParticleVariableCreate( Particles, Name, DOFs, Output, &
             Secondary, TYPE ) 
!------------------------------------------------------------------------------
    TYPE(Particle_t), POINTER :: Particles
    CHARACTER(LEN=*) :: Name
    INTEGER, OPTIONAL :: DOFs
    INTEGER, OPTIONAL :: TYPE
    LOGICAL, OPTIONAL :: Output
    LOGICAL, OPTIONAL :: Secondary
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Variables, Var
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER :: Dofs2
    LOGICAL :: stat
    INTEGER :: NoParticles
    TYPE(Mesh_t),POINTER :: Mesh
    TYPE(Solver_t),POINTER :: Solver
!------------------------------------------------------------------------------

    ! If already created don't do it again
    Var => VariableGet( Particles % Variables, Name )
    IF(ASSOCIATED(Var)) RETURN

    CALL Info('ParticleVariableCreate','Creating variable: '//TRIM(Name))

    NoParticles = Particles % MaxNumberOfParticles
    IF( NoParticles == 0 ) THEN
      CALL Warn('ParticleVariableCreate','No particles present!')
    END IF	

    IF( PRESENT( Dofs ) ) THEN
      Dofs2 = Dofs
    ELSE
      Dofs2 = 1
    END IF 

    NULLIFY( Values )
    ALLOCATE( Values( NoParticles * Dofs2 ) )
    Values = 0.0_dp

    Solver => CurrentModel % Solver
    Mesh => CurrentModel % Solver % Mesh

    CALL VariableAdd( Particles % Variables,Mesh,Solver,Name,DOFs2,Values,&
         Output=Output, Secondary=Secondary, TYPE=TYPE )

!------------------------------------------------------------------------------
  END SUBROUTINE ParticleVariableCreate
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Creates a variable related to the particles. The normal variabletype without
!> the permutation vector is used. 
!------------------------------------------------------------------------------
  SUBROUTINE ParticleVariableInitialize( Particles, Mesh, ToName, FromName )
!------------------------------------------------------------------------------
    TYPE(Particle_t), POINTER :: Particles
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: ToName
    CHARACTER(LEN=*), OPTIONAL :: FromName
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: ToVar, FromVar
    INTEGER :: Dofs,i,j,k,dim
    LOGICAL :: stat, UseHandle
    INTEGER :: NoParticles
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: val, vals(3), Coord(3), Basis(27), detJ
    TYPE(ValueHandle_t), SAVE :: InitField_h
!------------------------------------------------------------------------------

    NoParticles = Particles % NumberOfParticles
    

    k = COUNT( Particles % ElementIndex(1:NoParticles) == 0 )
    k = ParallelReduction(k)

    IF( k > 0 ) THEN
      CALL Info('ParticleVariableInitialize','We need owner elements, going to find them!',Level=10)
      CALL LocateParticles( Particles ) 
    END IF
        
    ToVar => VariableGet( Particles % Variables, TRIM(ToName) )
    IF(.NOT. ASSOCIATED(ToVar)) RETURN
    dofs = ToVar % Dofs
    
    UseHandle = .FALSE.
    IF( PRESENT(FromName)) THEN
      FromVar => VariableGet( Mesh % Variables, FromName )
    ELSE
      FromVar => VariableGet( Mesh % Variables, ToName )
    END IF
    IF(ASSOCIATED(FromVar)) THEN
      IF(dofs /= FromVar % Dofs) THEN
        CALL Fatal('ParticleVariableInitialize','To and From vars have different dofs!')
      END IF
      CALL Info('ParticleVariableInitialize','Initializing variable from variable in mesh: '//TRIM(ToName),Level=6)
    ELSE
      IF(dofs /= 1 ) THEN
        CALL Fatal('ParticleVariableInitialize','Initialization for scalars so far!')
      END IF
      CALL ListInitElementKeyword( InitField_h,'Initial Condition',ToName)
      IF( InitField_h % NotPresentAnywhere ) RETURN
      CALL Info('ParticleVariableInitialize','Initializing variable from initial section: '//TRIM(ToName),Level=6)
      UseHandle = .TRUE.
    END IF
      
    dim = Particles % dim
    
    DO i = 1, Particles % NumberOfParticles
      j = GetParticleElement( Particles, i )
      IF( j == 0 ) CYCLE

      Element => Mesh % Elements( j )
      CurrentModel % CurrentElement => Element
      Coord = 0._dp
      Coord(1:dim) = Particles % Coordinate(i, 1:dim) 

      stat = ParticleElementInfo( Element, Coord, DetJ, Basis )
      IF(.NOT. stat) CYCLE
      
      IF( UseHandle ) THEN
        val = ListGetElementReal( InitField_h, Basis, Element, Stat )
        ToVar % Values( i ) = val
      ELSE
        IF( ToVar % dofs == 1 ) THEN
          CALL GetScalarFieldInMesh(FromVar, Element, Basis, val ) 
          ToVar % Values( i ) =  val 
        ELSE
          CALL GetVectorFieldInMesh(FromVar, Element, Basis, vals ) 
          DO j=1,ToVar % dofs             
            ToVar % Values( dofs*(i-1)+j ) = vals(j)     
          END DO
        END IF
      END IF
    END DO

    IF( InfoActive(20) ) THEN
      PRINT *,'ParticleVariableInitialize',TRIM(ToVar % Name), MINVAL(ToVar % Values),MAXVAL(ToVar % Values)
    END IF
      
  END SUBROUTINE ParticleVariableInitialize
    
  
!------------------------------------------------------------------------------
!>  Given a variable name, get a handle to it.
!------------------------------------------------------------------------------
  FUNCTION ParticleVariableGet( Particles, Name ) RESULT ( Var )
    
    TYPE(Particle_t), POINTER :: Particles
    CHARACTER(LEN=*) :: Name
    TYPE(Variable_t), POINTER :: Var
!------------------------------------------------------------------------------
    Var => VariableGet( Particles % Variables, Name, ThisOnly = .TRUE. )

  END FUNCTION ParticleVariableGet
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Given a variable name, get a handle to its values.
!------------------------------------------------------------------------------
  FUNCTION ParticleVariableValues( Particles, Name ) RESULT ( Values )
    
    TYPE(Particle_t), POINTER :: Particles
    CHARACTER(LEN=*) :: Name
    REAL(KIND=dp), POINTER :: Values(:)
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var
    Var => VariableGet( Particles % Variables, Name )
    Values => Var % Values

  END FUNCTION ParticleVariableValues
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Resize the existing variables and optionally repack the old variables 
!> so that they lost ones are eliminated. This is controlled by the Perm
!> vector that includes only the indexes of the particles to be saved.
!------------------------------------------------------------------------------
    SUBROUTINE ParticleVariablesResize( Particles, Prevsize, Newsize, Perm) 
!------------------------------------------------------------------------------
      TYPE(Particle_t), POINTER :: Particles
      INTEGER :: Newsize, Prevsize
      INTEGER, OPTIONAL :: Perm(:)
!------------------------------------------------------------------------------
      TYPE(Variable_t), POINTER :: Var
      REAL(KIND=dp), POINTER :: Values(:)
      INTEGER :: i,j,k,OldSize
!------------------------------------------------------------------------------
      
      oldsize = prevsize
      IF( PRESENT( Perm ) ) THEN
        oldsize = SIZE(Perm)
        DO i=1,SIZE(Perm)
          IF( Perm(i) == 0 ) THEN
            oldsize = i-1
            EXIT
          END IF
        END DO
        IF( oldsize < 1 ) oldsize = 1
        CALL Info('ParticleVariablesResize','Using compact size of: '//I2S(oldsize),Level=12)
      END IF

      Var => Particles % Variables
      DO WHILE( ASSOCIATED(Var) )
        k = Var % NameLen
        IF( k > 0 ) THEN
          IF( Var % Dofs > 1 ) THEN
            CALL Fatal('ParticleVariableResize','Implement size increase for vectors!')
          END IF
          Values => Var % Values
          IF( SIZE( Var % Values ) < newsize ) THEN            
            CALL Info('ParticleVariableResize','Increasing size of variable: '// &
                Var % name(1:k),Level=12)
            ALLOCATE( Var % Values(newsize) )
            IF( PRESENT( Perm ) ) THEN
              Var % Values(1:oldsize) = Values(Perm(1:oldsize))
            ELSE
              OldSize = SIZE( Values ) 
              Var % Values(1:oldsize) = Values(1:oldsize)
            END IF
            Var % Values(oldsize+1:newsize) = 0.0_dp
            DEALLOCATE( Values ) 
          ELSE IF( PRESENT( Perm ) ) THEN
            CALL Info('ParticleVariableResize','Reorder dofs in variable: '// &
                Var % name(1:k),Level=15)
            Var % Values(1:oldsize) = Values(Perm(1:oldsize))
            Var % Values(oldsize+1:newsize) = 0.0_dp
          END IF
        END IF
        Var => Var % Next
      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE ParticleVariablesResize
!------------------------------------------------------------------------------


  
  !------------------------------------------------------------------------
  !> Write particles to an external file in simple ascii format matrix.
  !-------------------------------------------------------------------------
  SUBROUTINE ParticleOutputTable( Particles ) 
    
    TYPE(Particle_t), POINTER :: Particles
    
    TYPE(Variable_t), POINTER :: TimeVar
    TYPE(ValueList_t), POINTER :: Params 
    CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix, FileName
    LOGICAL :: Found, NumberFilesByParticles, NumberFilesBySteps, ParticleMode
    REAL(KIND=dp), POINTER :: Coord(:,:), Velo(:,:), Dist(:)
    REAL(KIND=dp) :: time
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Status(:)
    INTEGER :: i,j,n,dofs,Vari,Rank,dim, NoParticles, MinSaveStatus, MaxSaveStatus, WritePE
    INTEGER :: VisitedTimes = 0
    INTEGER, POINTER :: Indexes(:),Perm(:)
    LOGICAL :: GotTimeVar, GotDistVar
    TYPE(Variable_t), POINTER :: PartTimeVar, PartDistVar    
    REAL(KIND=dp), POINTER :: Basis(:)
    TYPE(Nodes_t) :: Nodes      
    INTEGER, PARAMETER :: TableUnit = 10


    SAVE :: VisitedTimes, Params, FilePrefix, NumberFilesByParticles, NumberFilesBySteps, &
        MinSaveStatus, MaxSaveStatus, TimeVar, Basis, Nodes

    CALL Info('ParticleOutputTable','Saving particle data into simple ascii table',Level=8)
    
    VisitedTimes = VisitedTimes + 1
    
    Mesh => GetMesh()
    dim = Particles % dim
    
    Coord => Particles % Coordinate
    Velo => Particles % Velocity 
    Status => Particles % Status
    
    PartDistVar => ParticleVariableGet( Particles,'particle distance')
    GotDistVar = ASSOCIATED( PartDistVar )

    PartTimeVar => ParticleVariableGet( Particles,'particle time')
    GotTimeVar = ASSOCIATED( PartTimeVar )

    ParticleMode = .NOT. ASSOCIATED( Particles % UVW ) 
    IF(.NOT. ParticleMode ) THEN
      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )
    END IF

    WritePE = ParEnv % MyPe
    IF( ParEnv % PEs > 1 ) THEN
      WritePE = ParallelReduction( WritePE, 1 )
    END IF
    
    IF( VisitedTimes == 1 ) THEN
      Params => ListGetSolverParams()
      FilePrefix = ListGetString(Params,'Filename Prefix')
      IF( ParEnv % MyPe == WritePE ) THEN
        CALL WriteParticleFileNames(FilePrefix, dim)
      END IF
        
      NumberFilesByParticles = ListGetLogical( Params,'Filename Particle Numbering',Found) 
      NumberFilesBySteps = ListGetLogical( Params,'Filename Timestep Numbering',Found) 
      IF( NumberFilesByParticles .AND. NumberFilesBySteps ) THEN
        CALL Fatal('ParticleOutputTable','Files may be numbered either by steps or particles')
      END IF
      
      MinSaveStatus = ListGetInteger( Params,'Min Status for Saving',Found)
      IF(.NOT. Found ) MinSaveStatus = PARTICLE_INITIATED
      
      MaxSaveStatus = ListGetInteger( Params,'Max Status for Saving',Found)
      IF(.NOT. Found ) MaxSaveStatus = PARTICLE_LOST-1
      
      TimeVar => VariableGet( Mesh % Variables,'time')

      IF( ParEnv % PEs > 1 ) THEN
        WRITE( FilePrefix,'(A,A,I4.4)' ) TRIM(FilePrefix),'par',ParEnv % MyPe + 1
      END IF
    END IF
    
    time = TimeVar % Values(1)
    NoParticles = Particles % NumberOfParticles

    CALL Info('ParticleOutputTable','Saving at maximum '//I2S(NoParticles)//' particles',Level=6)
    
    IF( NumberFilesByParticles ) THEN
      DO i = 1, NoParticles
        CALL OpenParticleFile(FilePrefix, i)
        IF ( Particles % Status(i) > MaxSaveStatus .OR. &
             Particles % Status(i) < MinSaveStatus )  CYCLE
        CALL WriteParticleLine( dim, i ) 
        CALL CloseParticleFile()
      END DO
    ELSE
      IF( NumberFilesBySteps ) THEN
        CALL OpenParticleFile(FilePrefix, VisitedTimes )
      ELSE
        CALL OpenParticleFile(FilePrefix, 0 )        
      END IF
      DO i = 1, NoParticles
        IF ( Particles % Status(i) > MaxSaveStatus .OR. &
             Particles % Status(i) < MinSaveStatus )  CYCLE
        CALL WriteParticleLine( dim, i )        
      END DO
      CALL CloseParticleFile()
    END IF

    IF( .NOT. ParticleMode ) THEN
      DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z )
    END IF


  CONTAINS

    !------------------------------------------------------------------------
    !> Write the names file for user information. Remember to update this if 
    !> SaveParticleStep is modified.
    !-------------------------------------------------------------------------
    SUBROUTINE WriteParticleFileNames( Prefix, Dim ) 
      
      CHARACTER(*) :: Prefix
      INTEGER :: dim

      CHARACTER(:), ALLOCATABLE :: FileName
      INTEGER :: i,j,dofs
      TYPE(Variable_t), POINTER :: Solution
      LOGICAL :: ComponentVector, ThisOnly = .TRUE.
      CHARACTER(:), ALLOCATABLE :: Txt, FieldName


      FileName =  TRIM(FilePrefix)//'.dat.names'
      
      OPEN (TableUnit, FILE=FileName )
      
      WRITE( TableUnit, '(A)' ) 'Variables in file: '//TRIM(FilePrefix)//'*.dat'
      WRITE( TableUnit, '(A,I2)' ) 'Dimension of particle set is',dim
      i = 1
      WRITE( TableUnit, '(I2,A)' )  i,': time'

      IF( NumberFilesBySteps ) THEN
        WRITE( TableUnit,'(I2, A)' ) i+1, ': particle id'
        i = i + 1
      ELSE IF( NumberFilesByParticles ) THEN
        WRITE( TableUnit,'(I2, A)' ) i+1, ': visited time'
        i = i + 1
      ELSE     
        WRITE( TableUnit,'(I2, A)' ) i+1, ': visited time'  
        WRITE( TableUnit,'(I2, A)' ) i+2, ': particle id'
        i = i + 2
      END IF
      
      WRITE( TableUnit, '(I2,A)' )  i+1,': Coordinate_1'
      WRITE( TableUnit, '(I2,A)' )  i+2,': Coordinate_2'      
      IF(dim == 3) WRITE( TableUnit, '(I2,A)' )  i+3,': Coordinate_3'
      i = i + DIM

      IF( ParticleMode ) THEN
        WRITE( TableUnit, '(I2,A)' )  i+1,': Velocity_1'
        WRITE( TableUnit, '(I2,A)' )  i+2,': Velocity_2'
        IF(dim == 3) WRITE( TableUnit, '(I2,A)' )  i+3,': Velocity_3'
        i = i + DIM
        IF( GotDistVar ) THEN
          WRITE( TableUnit, '(I2,A)' )  i+1,': Particle Distance'
          i = i + 1
        END IF
        IF( GotTimeVar ) THEN
          WRITE( TableUnit, '(I2,A)' )  i+1,': Particle Time'
          i = i + 1
        END IF
      ELSE  
              
        DO Rank = 1,2

          DO Vari = 1, 99
            
            IF( Rank == 1 ) THEN
              Txt = 'Scalar Field '//I2S(Vari)
            ELSE
              Txt = 'Vector Field '//I2S(Vari)
            END IF
            
            FieldName = ListGetString( Params, Txt, Found )
            IF(.NOT. Found) EXIT
            
            Solution => VariableGet( Mesh % Variables, &
                TRIM(FieldName),ThisOnly )
            ComponentVector = .FALSE.
            IF( Rank == 2 ) THEN
              IF(.NOT. ASSOCIATED(Solution)) THEN
                Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 1', ThisOnly )
                IF( ASSOCIATED(Solution)) THEN 
                  ComponentVector = .TRUE.
                ELSE
                  CALL Warn('WriteParticleLine', 'Nonexistent variable: '//TRIM(FieldName))
                  CYCLE
                END IF
              END IF
            ELSE
              IF(.NOT. ASSOCIATED(Solution)) THEN
                CALL Warn('WriteParticleLine','Nonexistent variable: '//TRIM(FieldName))
                CYCLE
              END IF
            END IF
            
            IF( ASSOCIATED(Solution % EigenVectors)) THEN
              CALL Warn('WriteParticleLine','Do the eigen values')
            END IF
            
            dofs = Solution % DOFs
                    
            IF( ComponentVector ) THEN
              Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly )
              IF( ASSOCIATED(Solution)) THEN
                dofs = 2 
                Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly )
                IF( ASSOCIATED(Solution)) THEN
                  dofs = 3 
                END IF
              END IF
            END IF
            
            IF( dofs == 1 ) THEN
              WRITE( TableUnit, '(I2,A)' )  i+1,': '//TRIM(FieldName) 
              i = i + 1
            ELSE
              DO j=1,dofs                
               WRITE( TableUnit, '(I2,A)' )  i+j,': '//TRIM(FieldName)//'_'//I2S(j)  
             END DO
             i = i + dofs
           END IF

         END DO
       END DO
     END IF
     
     CLOSE( TableUnit )
      
   END SUBROUTINE WriteParticleFileNames



    !------------------------------------------------------------------------
    !> Open a numbered file for each particle. These must be separate since the 
    !> number of steps for each particle may vary greatly
    !-------------------------------------------------------------------------
    SUBROUTINE OpenParticleFile( Prefix, FileNo ) 
      
      CHARACTER(LEN=*) :: Prefix
      INTEGER :: FileNo
      LOGICAL, SAVE :: Visited = .FALSE.
      CHARACTER(:), ALLOCATABLE :: FileName
      !CHARACTER(MAX_NAME_LEN) :: FileName

      IF( FileNo == 0 ) THEN
        FileName = TRIM(FilePrefix)//'.dat'
        IF( .NOT. Visited ) THEN
          CALL Info( 'ParticleOutputTable', 'Saving particle data to file: '//FileName, Level=4 )
        END IF
      ELSE
        IF ( FileNo==1 .AND.  .NOT. Visited ) THEN
          CALL Info( 'ParticleOutputTable', 'Saving particle data to files: '// &
                             TRIM(FilePrefix)//'_*.dat', Level=4 )
        END IF
        FileName=TRIM(FilePrefix)//'_'//i2s(fileno)//'.dat'
      END IF
      
      IF( VisitedTimes == 1 .OR. NumberFilesBySteps ) THEN
        OPEN (TableUnit, FILE=FileName )
        WRITE( TableUnit, '(A)', ADVANCE='no' ) ' ' ! delete old contents
      ELSE
        OPEN (TableUnit, FILE=FileName,POSITION='APPEND' )
      END IF

      Visited = .TRUE.

    END SUBROUTINE OpenParticleFile
    
    
    !------------------------------------------------------------------------
    !> Save one line in the particle file
    !-------------------------------------------------------------------------
    SUBROUTINE WriteParticleLine( Dim, No )
      INTEGER :: Dim, No
      TYPE(Variable_t), POINTER :: Solution
      REAL(KIND=dp), POINTER :: Values(:)
      REAL(KIND=dp) :: u,v,w,val,detJ,r(3)
      LOGICAL :: stat, ThisOnly=.TRUE., ComponentVector
      CHARACTER(:), ALLOCATABLE :: Txt, FieldName
      TYPE(Element_t), POINTER :: Element

      INTEGER, ALLOCATABLE :: ElemInd(:)

      WRITE( TableUnit,'(ES12.4)', ADVANCE = 'NO' ) time
        
      IF( NumberFilesBySteps ) THEN
        WRITE( TableUnit,'(I9)', ADVANCE = 'NO' ) No
      ELSE IF( NumberFilesByParticles ) THEN
        WRITE( TableUnit,'(I9)', ADVANCE = 'NO' ) VisitedTimes
      ELSE       
        WRITE( TableUnit,'(2I9)', ADVANCE = 'NO' ) VisitedTimes, No
      END IF
      
      IF( ParticleMode ) THEN
       
        IF( dim == 3 ) THEN
          WRITE(TableUnit,'(6ES16.7E3)',ADVANCE='NO') Coord(No,1:3), Velo(No,1:3)
        ELSE
          WRITE(TableUnit,'(4ES16.7E3)',ADVANCE='NO') Coord(No,1:2), Velo(No,1:2)
        END IF

        IF( GotDistVar ) WRITE (TableUnit,'(ES12.4)',ADVANCE='NO') PartDistVar % Values(No)
        IF( GotTimeVar ) WRITE (TableUnit,'(ES12.4)',ADVANCE='NO') PartTimeVar % Values(No)
        
      ELSE

        Element => Mesh % Elements( Particles % ElementIndex(No) )            
        n = Element % TYPE % NumberOfNodes
        Indexes => Element % NodeIndexes
        
        Nodes % x(1:n) = Mesh % Nodes % x( Indexes ) 
        Nodes % y(1:n) = Mesh % Nodes % y( Indexes ) 
        Nodes % z(1:n) = Mesh % Nodes % z( Indexes ) 

        ALLOCATE(ElemInd(Mesh % MaxElementDOFs))

        u = Particles % uvw(No,1)
        v = Particles % uvw(No,2)
        IF( dim == 3 ) THEN
          w = Particles % uvw(No,3)
        ELSE
          w = 0.0_dp
        END IF
        
        stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis)
        
        r(1) = SUM( Basis(1:n) * Nodes % x(1:n) )
        r(2) = SUM( Basis(1:n) * Nodes % y(1:n) )
        r(3) = SUM( Basis(1:n) * Nodes % z(1:n) )

        WRITE( TableUnit,'(2ES16.7E3)', ADVANCE='no') r(1:2)             
        IF( dim == 3 ) THEN
          WRITE( TableUnit,'(ES16.7E3)', ADVANCE='no') r(3)
        END IF


        DO Rank = 1,2

          DO Vari = 1, 99
            
            IF( Rank == 1 ) THEN
              Txt = 'Scalar Field '//I2S(Vari)
            ELSE
              Txt = 'Vector Field '//I2S(Vari)
            END IF
            
            FieldName = ListGetString( Params, TRIM(Txt), Found )
            IF(.NOT. Found) EXIT
            
            Solution => VariableGet( Mesh % Variables, &
                TRIM(FieldName),ThisOnly )
            ComponentVector = .FALSE.
            IF( Rank == 2 ) THEN
              IF(.NOT. ASSOCIATED(Solution)) THEN
                Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 1', ThisOnly )
                IF( ASSOCIATED(Solution)) THEN 
                  ComponentVector = .TRUE.
                ELSE
                  CALL Warn('WriteParticleLine', 'Nonexistent variable: '//TRIM(FieldName))
                  CYCLE
                END IF
              END IF
            ELSE
              IF(.NOT. ASSOCIATED(Solution)) THEN
                CALL Warn('WriteParticleLine','Nonexistent variable: '//TRIM(FieldName))
                CYCLE
              END IF
            END IF
            
            IF( ASSOCIATED(Solution % EigenVectors)) THEN
              CALL Warn('WriteParticleLine','Do the eigen values')
            END IF
            
            Perm => Solution % Perm
            dofs = Solution % DOFs
            Values => Solution % Values         
            val = 0.0_dp

            IF( Solution % TYPE == Variable_on_nodes_on_elements ) THEN
              ElemInd(1:n) = Perm(Element % DGIndexes(1:n))
            ELSE
              ElemInd(1:n) = Perm(Indexes(1:n))
            END IF
            
            IF( ComponentVector ) THEN
              IF( ALL(ElemInd(1:n) > 0 ) ) THEN
                val = SUM( Basis(1:n) * Values(ElemInd(1:n)) )
              END IF
              WRITE( TableUnit,'(ES16.7E3)', ADVANCE='no') val             

              Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly )
              IF( ASSOCIATED(Solution)) THEN
                Values => Solution % Values
                IF( ALL( ElemInd(1:n) > 0 ) ) THEN
                  val = SUM( Basis(1:n) * Values(ElemInd(1:n)) )
                END IF
                WRITE( TableUnit,'(ES16.7E3)', ADVANCE='no') val             

                Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly )
                IF( ASSOCIATED(Solution)) THEN
                  Values => Solution % Values
                  IF( ALL( ElemInd(1:n) > 0 ) ) THEN
                    val = SUM( Basis(1:n) * Values(ElemInd(1:n)) )
                  END IF
                  WRITE( TableUnit,'(ES16.7E3)', ADVANCE='no') val             
                END IF
              END IF
            ELSE IF( Dofs == 1 ) THEN
              IF( ALL( ElemInd(1:n) > 0 ) ) THEN
                val = SUM( Basis(1:n) * Values(ElemInd(1:n)) )
              END IF
              WRITE( TableUnit,'(ES16.7E3)', ADVANCE='no') val             
            ELSE
              DO j = 1, Dofs
                IF( ALL( ElemInd(1:n) > 0 ) ) THEN
                  val = SUM( Basis(1:n) * Values(Dofs*(ElemInd(1:n)-1)+j) )
                END IF
                WRITE( TableUnit,'(ES16.7E3)', ADVANCE='no') val                          
              END DO
            END IF
          END DO
        END DO
        WRITE(TableUnit,'(A)') ' '
      END IF

    END SUBROUTINE WriteParticleLine
    
    
    !------------------------------------------------------------------------
    ! Close the particle file
    !-------------------------------------------------------------------------
    SUBROUTINE CloseParticleFile( )
      
      CLOSE( TableUnit )
      
    END SUBROUTINE CloseParticleFile
           
  END SUBROUTINE ParticleOutputTable



  !------------------------------------------------------------------------
  !> Write particles to an external file in Gmsh format as vector points.
  ! Subroutine contributed by Emilie Marchandise.
  !-------------------------------------------------------------------------
  SUBROUTINE ParticleOutputGmsh( Particles )
    
    TYPE(Particle_t), POINTER :: Particles
    
    TYPE(Variable_t), POINTER :: TimeVar
    TYPE(ValueList_t), POINTER :: Params 
    CHARACTER(:), ALLOCATABLE :: FilePrefix, FileNameGmsh, FileNameOut
    LOGICAL :: Found
    REAL(KIND=dp), POINTER :: Coord(:,:), Velo(:,:), Dist(:)
    REAL(KIND=dp), POINTER :: CoordInit(:,:)
    REAL(KIND=dp) :: time
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Status(:)
    INTEGER :: i,j,k, dim, NoParticles, nStep
    INTEGER :: VisitedTimes = 0
    INTEGER, POINTER :: TimeSteps(:)
    LOGICAL :: GotTimeVar, GotDistVar
    TYPE(Variable_t), POINTER :: PartTimeVar, PartDistVar    

    SAVE :: VisitedTimes, Params, FilePrefix, TimeVar, FileNameGmsh, CoordInit

    Params => ListGetSolverParams()
    FilePrefix = ListGetString(Params,'Filename Prefix')

    FileNameGmsh = FilePrefix // '.pos'
    
    Mesh => GetMesh()
    dim = Particles % dim

    Coord => Particles % Coordinate
    Velo => Particles % Velocity 
    Status => Particles % Status
    
    NoParticles = Particles % NumberOfParticles

    VisitedTimes = VisitedTimes + 1
    IF (VisitedTimes==1) THEN

       OPEN (10, FILE=FileNameGmsh )
       WRITE( 10, '(A)') 'View[0].VectorType=5; //for displacement type'
       WRITE( 10, '(A)') 'View[0].PointType=1; // for spheres'
       WRITE( 10, '(A)') 'View[0].PointSize=5; // for spheres'
       WRITE( 10, '(A)') 'View[0].IntervalsType = 1; //for iso-values interval'
       WRITE( 10, '(A)') 'View[0].NbIso = 1; //for one color'
       WRITE( 10, '(A)') 'View[0].ShowScale = 0; ' 
       CLOSE( 10 )

       ALLOCATE( CoordInit(NoParticles,dim) )
       CoordInit = Particles % Coordinate
  
       TimeVar => VariableGet( Mesh % Variables,'time')
    END IF

    time = TimeVar % Values(1)
    CALL Info( 'ParticleTracker', 'Saving particle paths to file: '//TRIM(FileNameGmsh), Level=4 )
    
    OPEN (10, FILE=FileNameGmsh, POSITION='APPEND' )   
    WRITE( 10, '(A)') 'View "particles" {'
    WRITE( 10, '(A)', ADVANCE='NO') 'TIME{'
    WRITE( 10, '(ES16.7E3)', ADVANCE='NO') time
    WRITE( 10, '(A)') '};'

    DO i = 1, NoParticles
        WRITE( 10, '(A)', ADVANCE='NO') 'VP('
        DO k=1,dim
           WRITE( 10, '(ES16.7E3)', ADVANCE='NO') CoordInit(i,k)
           IF(k < dim) WRITE( 10, '(A)', ADVANCE='NO') ','
        END DO
        IF (dim ==2)  WRITE( 10, '(A)', ADVANCE='NO') ', 0.0'
        WRITE( 10, '(A)', ADVANCE='NO') '){'
        DO k=1,dim
           WRITE( 10, '(ES16.7E3)', ADVANCE='NO') Coord(i,k)-CoordInit(i,k)
           IF(k < dim) WRITE( 10, '(A)', ADVANCE='NO') ','
        END DO
        IF (dim ==2)  WRITE( 10, '(A)', ADVANCE='NO') ', 0.0'
        WRITE( 10, '(A)') '};'
    END DO
    WRITE( 10, '(A)') '};'


    ! Save for last timestep, this is a conservative estimate assuming 
    ! savig after each timestep.
    !-----------------------------------------------------------------
    Timesteps => ListGetIntegerArray( CurrentModel % Simulation, &
        'Timestep Intervals', Found )
    nStep = SUM( TimeSteps )

    IF (VisitedTimes == nStep) THEN 
       WRITE( FileNameOut,'(A,A)') TRIM(FilePrefix),'_combined.pos";'
       WRITE( 10, '(A)') 'Combine TimeStepsByViewName;'
       WRITE( 10, '(A)', ADVANCE='NO') 'Save View [0] "'
       WRITE( 10, '(A)') FileNameOut
       WRITE( 10, '(A)') 'Printf("View[0].VectorType=5;'
       WRITE( 10, '(A)') 'View[0].PointType=1; // for spheres'
       WRITE( 10, '(A)') 'View[0].PointSize=5; // for spheres'
       WRITE( 10, '(A)') 'View[0].IntervalsType = 1; //for iso-values interval'
       WRITE( 10, '(A)') 'View[0].NbIso = 1; //for one color'
       WRITE( 10, '(A)', ADVANCE='NO') 'View[0].ShowScale = 0; ") >> "'
       WRITE( 10, '(A)') FileNameOut
    END IF

    CLOSE( 10 )


  END SUBROUTINE ParticleOutputGmsh


!------------------------------------------------------------------------------
!> Saves particles in unstructured XML VTK format (VTU) to an external file. 
!------------------------------------------------------------------------------
  SUBROUTINE ParticleOutputVtu( Particles )
!------------------------------------------------------------------------------

    USE DefUtils 
    USE MeshUtils
    USE ElementDescription
    USE AscBinOutputUtils
    
    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles  
    
    TYPE(ValueList_t),POINTER :: Params
    INTEGER, SAVE :: nTime = 0
    LOGICAL :: GotIt, Parallel, FixedMeshend,SinglePrec
    
    CHARACTER(:), ALLOCATABLE :: FilePrefix, BaseFile, OutputDirectory
    CHARACTER(MAX_NAME_LEN) :: VtuFile, PvtuFile
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i, j, k, Partitions, Part, ExtCount, FileindexOffSet, &
        Status, MinSaveStatus, MaxSaveStatus, PrecBits, PrecSize, IntSize, &
        iTime, SaveGroup 
    REAL(KIND=dp) :: SaveNodeFraction, LocalVal(3)
    LOGICAL :: BinaryOutput,AsciiOutput,Found,Visited = .FALSE.,SaveFields
    REAL(KIND=dp) :: DoubleWrk
    REAL :: SingleWrk

    TYPE(Solver_t), POINTER :: pSolver
    INTEGER :: NumberOfNodes, ParallelNodes, Dim
    
    SAVE :: MinSaveStatus, MaxSaveStatus, FilePrefix
    
    Params => ListGetSolverParams()
    Mesh => GetMesh()
    pSolver => GetSolver()
    
    ExtCount = ListGetInteger( Params,'Output Count',GotIt)
    IF( GotIt ) THEN
      nTime = ExtCount
    ELSE
      nTime = nTime + 1
    END IF
    FileIndexOffset = ListGetInteger( Params,'Fileindex offset',GotIt)
    iTime = nTime + FileIndexOffset

    IF ( nTime == 1 ) THEN
      FilePrefix = ListGetString( Params,'Filename Prefix')
      CALL Info('ParticleOutputVtu','Saving in VTK XML unstructured format to file: ' &
	//TRIM(FilePrefix)//'.vtu')
      
      MinSaveStatus = ListGetInteger( Params,'Min Status for Saving',Found)
      IF(.NOT. Found ) MinSaveStatus = PARTICLE_INITIATED
      
      MaxSaveStatus = ListGetInteger( Params,'Max Status for Saving',Found)
      IF(.NOT. Found ) MaxSaveStatus = PARTICLE_LOST-1      
    END IF
    
    BinaryOutput = GetLogical( Params,'Binary Output',GotIt)
    IF( GotIt ) THEN
      AsciiOutput = .NOT. BinaryOutput
    ELSE
      AsciiOutput = GetLogical( Params,'Ascii Output',GotIt)
      BinaryOutput = .NOT. AsciiOutput
    END IF
    
    SaveFields = GetLogical( Params,'Save Fields',GotIt)
    IF(.NOT. GotIt) SaveFields = .TRUE.

    SinglePrec = GetLogical( Params,'Single Precision',GotIt) 
    IF( SinglePrec ) THEN
      CALL Info('ParticleOutputVtu','Using single precision arithmetics in output!',Level=7)
    END IF
    
    IF( SinglePrec ) THEN
      PrecBits = 32
      PrecSize = KIND( SingleWrk ) 
    ELSE
      PrecBits = 64
      PrecSize = KIND( DoubleWrk ) 
    END IF
    IntSize = KIND(i)
    
    Partitions = ParEnv % PEs
    Part = ParEnv % MyPE
    Parallel = (Partitions > 1) .OR. GetLogical(Params,'Enforce Parallel format',GotIt)
    
    Dim = Particles % dim

    SaveGroup = ListGetInteger( Params,'Particle Save Group',GotIt)

    NumberOfNodes = 0
    DO i=1,Particles % NumberOfParticles
      IF ( Particles % Status(i) > MaxSaveStatus .OR. &
          Particles % Status(i) < MinSaveStatus )  CYCLE
      IF(SaveGroup > 0) THEN
        IF( Particles % Group(i) /= SaveGroup ) CYCLE 
      END IF    
      NumberOfNodes = NumberOfNodes + 1
    END DO

    SaveNodeFraction = ListGetCReal( Params,'Particle Save Fraction',GotIt)
    IF(GotIt) THEN
      NumberOfNodes = NINT( SaveNodeFraction * NumberOfNodes )
    ELSE
      i = ListGetInteger( Params,'Particle Save Number',GotIt)
      IF( GotIt ) THEN
        NumberOfNodes = MIN(i,NumberOfNodes)
      END IF
    END IF
    
    BaseFile = FilePrefix
    CALL SolverOutputDirectory( pSolver, BaseFile, OutputDirectory, UseMeshDir = .TRUE.  )
    BaseFile = TRIM(OutputDirectory)// '/' //TRIM(BaseFile)    

    IF(Parallel .AND. Part == 0) THEN
      IF( iTime < 10000 ) THEN
        WRITE( PvtuFile,'(A,I4.4,".pvtu")' ) TRIM(BaseFile),iTime
      ELSE
        WRITE( PvtuFile,'(A,I0,".pvtu")' ) TRIM(BaseFile),iTime
      END IF
      CALL WritePvtuFile( PvtuFile )
    END IF
    
    IF ( Parallel ) THEN
      IF( iTime < 10000 ) THEN
        WRITE( VtuFile,'(A,I4.4,A,I4.4,".vtu")' ) TRIM(BaseFile),Part+1,"par",&
            iTime
      ELSE
        WRITE( VtuFile,'(A,I4.4,A,I0,".vtu")' ) TRIM(BaseFile),Part+1,"par",&
            iTime
      END IF
    ELSE
      IF( iTime < 10000 ) THEN
        WRITE( VtuFile,'(A,I4.4,".vtu")' ) TRIM(BaseFile),iTime
      ELSE
        WRITE( VtuFile,'(A,I0,".vtu")' ) TRIM(BaseFile),iTime
      END IF
    END IF

    CALL Info('ParticleOutputVtu','Saving particles to file: '//TRIM(VtuFile),Level=8)
    CALL WriteVtuFile( VtuFile )
    

  CONTAINS

    
    SUBROUTINE WriteVtuFile( VtuFile )
      CHARACTER(LEN=*), INTENT(IN) :: VtuFile
      INTEGER, PARAMETER :: VtuUnit = 58
      TYPE(Variable_t), POINTER :: Var, Solution
      INTEGER :: i,j,k,dofs,Rank,cumn,n,vari,sdofs,IsVector,Offset,PartDim

      CHARACTER(:), ALLOCATABLE :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
          FieldName2, BaseString

      CHARACTER(1024) :: OutStr

      CHARACTER :: lf
      LOGICAL :: ScalarsExist, VectorsExist, Found, ParticleMode, ComponentVector, &
          ComplementExists, ThisOnly, Stat
      LOGICAL :: WriteData, WriteXML, Buffered, IsDG, IsInteger
      INTEGER, POINTER :: Perm(:), Perm2(:), Indexes(:)
      INTEGER, ALLOCATABLE :: ElemInd(:),ElemInd2(:)
      REAL(KIND=dp), POINTER :: Values(:),Values2(:),&
          Values3(:),VecValues(:,:),Basis(:)
      INTEGER, POINTER :: Ivalues(:)
      REAL(KIND=dp) :: x,y,z,u,v,w,DetJ,val
      TYPE(Nodes_t) :: Nodes      
      TYPE(Element_t), POINTER :: Element
      TYPE(Variable_t), POINTER :: ParticleVar

      ! Initialize the auxiliary module for buffered writing
      !--------------------------------------------------------------
      CALL AscBinWriteInit( AsciiOutput, SinglePrec, VtuUnit, NumberOfNodes )

      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

      n = Mesh % MaxElementDOFS
      ALLOCATE( ElemInd(n), ElemInd2(n) )

      ThisOnly = .TRUE.

      ParticleMode = .NOT. ASSOCIATED( Particles % UVW ) 

      ! Linefeed character
      !-----------------------------------
      lf = CHAR(10)
      dim = 3

      PartDim = Particles % Dim
      
      WriteXML = .TRUE.
      WriteData = AsciiOutput
      Params => ListGetSolverParams()
      Buffered = .TRUE.
      
      ! This is a hack to ensure that the streamed saving will cover the whole file
      !----------------------------------------------------------------------------
      IF(.TRUE.) THEN
        OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'formatted', STATUS='unknown' )
        WRITE( VtuUnit,'(A)') ' '
        CLOSE( VtuUnit ) 
      END IF
      
      ! This format works both for ascii and binary output
      !-------------------------------------------------------------------------
      OPEN( UNIT=VtuUnit, FILE=VtuFile, FORM = 'unformatted', ACCESS = 'stream', STATUS='unknown' )
      
      WRITE( OutStr,'(A)') '<?xml version="1.0"?>'//lf
      CALL AscBinStrWrite( OutStr ) 

      IF ( LittleEndian() ) THEN
        OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
      ELSE
        OutStr = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'//lf
      END IF
      CALL AscBinStrWrite( OutStr )
      WRITE( OutStr,'(A)') '  <UnstructuredGrid>'//lf
      CALL AscBinStrWrite( OutStr )
      WRITE( OutStr,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="',NumberOfNodes,&
          '" NumberOfCells="',NumberOfNodes,'">'//lf
      CALL AscBinStrWrite( OutStr )
    
      !---------------------------------------------------------------------
      ! Header information for nodewise data
      !---------------------------------------------------------------------
      ScalarFieldName = ListGetString( Params,'Scalar Field 1',ScalarsExist)
      VectorFieldName = ListGetString( Params,'Vector Field 1',VectorsExist)
      IF( .NOT. ( ScalarsExist .OR. VectorsExist) ) THEN
        CALL Warn('WriteVtuFile','Are there really no scalars or vectors?')
      END IF

      WRITE( OutStr,'(A)') '      <PointData>'//lf
      CALL AscBinStrWrite( OutStr )
      
      !---------------------------------------------------------------------
      ! do the scalars & vectors
      !--------------------------------- -----------------------------------
100   Offset = 0
      IsInteger = .FALSE.

      IF( SaveFields ) THEN
        
        DO IsVector = 0, 1

          DO Vari = 1, 999

            IF( IsVector == 0 ) THEN
              BaseString = 'Scalar Field'
            ELSE
              BaseString = 'Vector Field'
            END IF

            Txt = TRIM(BaseString)//' '//I2S(Vari)
            FieldName = ListGetString( Params, TRIM(Txt), Found )
            IF(.NOT. Found) EXIT

            !---------------------------------------------------------------------
            ! Find the variable with the given name in the normal manner 
            !---------------------------------------------------------------------
            IF( .NOT. ParticleMode ) THEN
              Solution => VariableGet( Mesh % Variables,TRIM(FieldName),ThisOnly )
              IF( ASSOCIATED( Solution ) ) THEN
                IF( IsVector == 1 .AND. Solution % Dofs <= 1 ) NULLIFY( Solution ) 
              END IF

              ComponentVector = .FALSE.
              IF( IsVector == 1 ) THEN
                IF(.NOT. ASSOCIATED(Solution)) THEN
                  Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 1', ThisOnly )
                  ComponentVector = ASSOCIATED( Solution ) 
                END IF
              END IF
              IF( .NOT. ASSOCIATED( Solution ) ) THEN
                CALL Warn('WriteVtuXMLFile','Nonexistent variable: '//TRIM(FieldName))
                CYCLE
              END IF

              IF( ASSOCIATED(Solution % EigenVectors)) THEN
                CALL Warn('WriteVtuXMLFile','Do the eigen values')
              END IF

              Perm => Solution % Perm
              dofs = Solution % DOFs
              Values => Solution % Values

              !---------------------------------------------------------------------
              ! Some vectors are defined by a set of components (either 2 or 3)
              !---------------------------------------------------------------------
              IF( ComponentVector ) THEN
                Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly )
                IF( ASSOCIATED(Solution)) THEN
                  Values2 => Solution % Values
                  dofs = 2
                END IF
                Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly )
                IF( ASSOCIATED(Solution)) THEN
                  Values3 => Solution % Values
                  dofs = 3
                END IF
              END IF

              !---------------------------------------------------------------------
              ! There may be special complementary variables such as 
              ! displacement & mesh update 
              !---------------------------------------------------------------------
              ComplementExists = .FALSE.
              Txt = TRIM(BaseString)//' '//I2S(Vari)//' Complement'

              FieldName2 = ListGetString( Params, TRIM(Txt), Found )
              IF( Found ) THEN
                Solution => VariableGet( Mesh % Variables, &
                    TRIM(FieldName2), ThisOnly )
                IF( ASSOCIATED(Solution)) THEN 
                  Values2 => Solution % Values
                  Perm2 => Solution % Perm 
                  ComplementExists = .TRUE.
                ELSE
                  CALL Warn('WriteVTUFile','Complement does not exist:'//TRIM(FieldName2))
                END IF
              END IF
            END IF

            !---------------------------------------------------------------------
            ! Get the values assuming particle mode
            !---------------------------------------------------------------------
            IF( ParticleMode ) THEN
              IsInteger = .FALSE.
              IF( IsVector == 1) THEN
                dofs = PartDim
                IF( FieldName == 'velocity' ) THEN
                  VecValues => Particles % Velocity
                ELSE IF( FieldName == 'force') THEN
                  VecValues => Particles % Force 
                ELSE
                  CALL Warn('WriteVtuXMLFile', 'Nonexistent variable: '//TRIM(FieldName))
                  CYCLE
                END IF
              ELSE
                dofs = 1
                IF( FieldName == 'particle distance') THEN
                  ParticleVar => ParticleVariableGet( Particles,'particle distance' )
                  IF( .NOT. ASSOCIATED( ParticleVar ) ) THEN
                    CALL Fatal('WriteVTUFile','> Particle Distance < does not exist!')
                  END IF
                  Values => ParticleVar % Values
                ELSE IF( FieldName == 'particle time') THEN
                  ParticleVar => ParticleVariableGet( Particles,'particle time' )
                  IF( .NOT. ASSOCIATED( ParticleVar ) ) THEN
                    CALL Fatal('WriteVTUFile','> Particle Time < does not exist!')
                  END IF
                  Values => ParticleVar % Values
                ELSE IF( FieldName == 'particle dt') THEN
                  ParticleVar => ParticleVariableGet( Particles,'particle dt' )
                  IF( .NOT. ASSOCIATED( ParticleVar ) ) THEN
                    CALL Fatal('WriteVTUFile','> Particle Dt < does not exist!')
                  END IF
                  Values => ParticleVar % Values
                ELSE IF( FieldName == 'particle group') THEN
                  IsInteger = .TRUE.
                  IValues => Particles % Group 
                  IF( .NOT. ASSOCIATED( IValues) ) THEN
                    CALL Fatal('WriteVTUFile','> Particle Group < does not exist!')
                  END IF
                ELSE
                  CALL Warn('WriteVtuXMLFile','Nonexistent variable: '//TRIM(FieldName))
                  CYCLE
                END IF
              END IF
            END IF

            !---------------------------------------------------------------------
            ! Finally save the field values for scalars
            !---------------------------------------------------------------------
            IF( dofs == 1 ) THEN
              sdofs = 1
            ELSE
              sdofs = MAX(dofs,3)
            END IF


            IF( WriteXML ) THEN
              WRITE( OutStr,'(A,I0,A)') '        <DataArray type="Float',PrecBits,'" Name="'//TRIM(FieldName)
              CALL AscBinStrWrite( OutStr )

              WRITE( OutStr,'(A,I0,A)') '" NumberOfComponents="',sdofs,'"'          
              CALL AscBinStrWrite( OutStr ) 

              IF( AsciiOutput ) THEN
                WRITE( OutStr,'(A)') ' format="ascii">'//lf
                CALL AscBinStrWrite( OutStr ) 
              ELSE
                WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
                CALL AscBinStrWrite( OutStr ) 
              END IF
            END IF


            IF( BinaryOutput ) THEN
              k = NumberOfNodes * PrecSize * sdofs
              Offset = Offset + IntSize + k
            END IF

            j = 0
            IF( WriteData ) THEN
              IF( BinaryOutput ) WRITE( VtuUnit ) k

              DO i = 1, Particles % NumberOfParticles

                IF ( Particles % Status(i) > MaxSaveStatus .OR. &
                    Particles % Status(i) < MinSaveStatus )  CYCLE
                IF( SaveGroup > 0 ) THEN
                  IF( Particles % Group(i) /= SaveGroup ) CYCLE
                END IF
                
                
                j = j + 1

                LocalVal = 0.0_dp

                IF( ParticleMode ) THEN
                  IF( IsVector == 1) THEN
                    dofs = dim
                    LocalVal(1:dofs) = VecValues(i,1:dim)
                  ELSE
                    dofs = 1
                    IF( IsInteger ) THEN
                      LocalVal(1) = 1.0_dp * IValues(i)
                    ELSE
                      LocalVal(1) = Values(i)
                    END IF
                  END IF
                ELSE
                  Element => Mesh % Elements( Particles % ElementIndex(i) )            
                  n = Element % TYPE % NumberOfNodes
                  Indexes => Element % NodeIndexes

                  Nodes % x(1:n) = Mesh % Nodes % x( Indexes ) 
                  Nodes % y(1:n) = Mesh % Nodes % y( Indexes ) 
                  Nodes % z(1:n) = Mesh % Nodes % z( Indexes ) 

                  u = Particles % uvw(i,1)
                  v = Particles % uvw(i,2)
                  IF( dim == 3 ) THEN
                    w = Particles % uvw(i,3)
                  ELSE
                    w = 0.0_dp
                  END IF

                  stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis)

                  IsDG = .FALSE.
                  IF( ASSOCIATED( Solution ) ) THEN
                    IsDG = ( Solution % TYPE == Variable_on_nodes_on_elements )
                  END IF
                  
                  IF( IsDG ) THEN
                    ElemInd(1:n) = Perm( Element % DGIndexes(1:n) )
                    IF( ComplementExists ) THEN
                      ElemInd2(1:n) = Perm2( Element % DGIndexes(1:n) )
                    END IF
                  ELSE 
                    ElemInd(1:n) = Perm( Indexes(1:n) )
                    IF( ComplementExists ) THEN
                      ElemInd2(1:n) = Perm2( Indexes(1:n) )
                    END IF
                  END IF

                  DO k=1,sdofs                  
                    val = 0.0_dp
                    IF( k <= dofs ) THEN
                      IF( ComponentVector ) THEN
                        IF( ALL( Perm( Indexes ) > 0 ) ) THEN
                          IF( k == 1 ) THEN
                            LocalVal(1) = SUM( Basis(1:n) * Values(ElemInd(1:n)) )
                          ELSE IF( k == 2 ) THEN
                            LocalVal(2) = SUM( Basis(1:n) * Values2(ElemInd(1:n)) )
                          ELSE
                            LocalVal(3) = SUM( Basis(1:n) * Values3(ElemInd(1:n)) )
                          END IF
                        END IF
                      ELSE 
                        IF( ALL( ElemInd(1:n) > 0 ) ) THEN
                          LocalVal(k) = SUM( Basis(1:n) * Values(dofs*(ElemInd(1:n)-1)+k) )
                        ELSE IF ( ComplementExists ) THEN
                          IF( ALL( ElemInd2(1:n) > 0 ) ) THEN
                            LocalVal(k) = SUM( Basis(1:n) * Values2(dofs*(ElemInd2(1:n)-1)+k) )
                          END IF
                        END IF
                      END IF
                    END IF
                  END DO
                END IF
                
                DO k=1,sdofs
                  CALL AscBinRealWrite( LocalVal(k) )
                END DO

                IF( j == NumberOfNodes ) EXIT
              END DO
              CALL AscBinRealWrite( 0.0_dp, .TRUE.)
            END IF

            IF( AsciiOutput ) THEN
              WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
              CALL AscBinStrWrite( OutStr ) 
            END IF
          END DO
        END DO
      END IF

      

      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '      </PointData>'//lf
        CALL AscBinStrWrite( OutStr ) 

        WRITE( OutStr,'(A)') '      <CellData>'//lf
        CALL AscBinStrWrite( OutStr ) 
        WRITE( OutStr,'(A)') '      </CellData>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
     

      ! Coordinates of each point
      !-------------------------------------
      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '      <Points>'//lf
        CALL AscBinStrWrite( OutStr ) 
        
        WRITE( OutStr,'(A,I0,A,I0,A)') '        <DataArray type="Float',PrecBits,'" NumberOfComponents="',dim,'"'
        CALL AscBinStrWrite( OutStr )       
        
        IF( AsciiOutput ) THEN
          WRITE( OutStr,'(A)') ' format="ascii">'//lf
          CALL AscBinStrWrite( OutStr ) 
        ELSE
          WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF


      IF( BinaryOutput ) THEN
        k = dim * NumberOfNodes * PrecSize
        Offset = Offset + IntSize + k
      END IF
      
      IF( WriteData ) THEN        
        IF( BinaryOutput ) WRITE( VtuUnit ) k

        LocalVal = 0.0_dp
        j = 0
        DO i = 1, Particles % NumberOfParticles
          
          IF ( Particles % Status(i) > MaxSaveStatus .OR. &
              Particles % Status(i) < MinSaveStatus )  CYCLE
          IF( SaveGroup > 0 ) THEN
            IF( Particles % Group(i) /= SaveGroup ) CYCLE
          END IF
          
          j = j + 1
          
          IF( ParticleMode ) THEN
            DO k=1,PartDim
              LocalVal(k) = Particles % Coordinate(i,k)
            END DO
          ELSE
            Element => Mesh % Elements( Particles % ElementIndex(i) )
            Indexes => Element % NodeIndexes
            n = Element % TYPE % NumberOfNodes
            
            Nodes % x(1:n) = Mesh % Nodes % x( Indexes ) 
            Nodes % y(1:n) = Mesh % Nodes % y( Indexes ) 
            Nodes % z(1:n) = Mesh % Nodes % z( Indexes ) 
            
            u = Particles % uvw(i,1)
            v = Particles % uvw(i,2)
            IF( dim == 3 ) THEN
              w = Particles % uvw(i,3)
            ELSE
              w = 0.0_dp
            END IF
            
            stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis)
            
            LocalVal(1) = SUM( Basis(1:n) * Nodes % x(1:n) )
            LocalVal(2) = SUM( Basis(1:n) * Nodes % y(1:n) )
            LocalVal(3)  = SUM( Basis(1:n) * Nodes % z(1:n) )
          END IF

          CALL AscBinRealWrite( LocalVal(1) )
          CALL AscBinRealWrite( LocalVal(2) )
          CALL AscBinRealWrite( LocalVal(3) )

          IF( j == NumberOfNodes ) EXIT
        END DO
        
        CALL AscBinRealWrite( 0.0_dp, .TRUE.)
      END IF

      IF( AsciiOutput ) THEN   
        WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '      </Points>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF


      ! Write out the mesh
      !-------------------------------------
      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '      <Cells>'//lf
        CALL AscBinStrWrite( OutStr ) 
        
        WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="connectivity"'
        CALL AscBinStrWrite( OutStr ) 
        
        IF( AsciiOutput ) THEN
          WRITE( OutStr,'(A)') ' format="ascii">'//lf
          CALL AscBinStrWrite( OutStr ) 
        ELSE
          WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF
      
      IF( BinaryOutput ) THEN
        ! The offset needs to be summed over all nodes
        k = NumberOfNodes * IntSize
        Offset = Offset + k + IntSize
      END IF
            
      IF( WriteData ) THEN
        IF( BinaryOutput ) WRITE( VtuUnit ) k        
        DO i = 1, NumberOfNodes
          CALL AscBinIntegerWrite( i - 1)
        END DO
        CALL AscBinIntegerWrite( 0, .TRUE. ) 
      END IF

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF

      ! Offsets for element indexes 
      !-------------------------------------------------------------------
      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="offsets"'
        CALL AscBinStrWrite( OutStr ) 
        
        IF( AsciiOutput ) THEN
          WRITE( OutStr,'(A)') ' format="ascii">'//lf
          CALL AscBinStrWrite( OutStr ) 
        ELSE
          WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF
      
      IF( BinaryOutput ) THEN
        k = NumberOfNodes * IntSize
        Offset = Offset + IntSize + k
      END IF
      
      IF( WriteData ) THEN
        IF( BinaryOutput ) WRITE( VtuUnit ) k         
        DO i = 1, NumberOfNodes
          CALL AscBinIntegerWrite( i )
        END DO
        CALL AscBinIntegerWrite( 0, .TRUE.)        
      END IF

      IF( AsciiOutput ) THEN   
        WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
      IF( WriteXML ) THEN
        WRITE( OutStr,'(A)') '        <DataArray type="Int32" Name="types"'
        CALL AscBinStrWrite( OutStr ) 
        
        IF( AsciiOutput ) THEN
          WRITE( OutStr,'(A)') ' FORMAT="ascii">'//lf
          CALL AscBinStrWrite( OutStr ) 
        ELSE
          WRITE( OutStr,'(A,I0,A)') ' format="appended" offset="',Offset,'"/>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF

      IF( BinaryOutput ) THEN
        k = NumberOfNodes * IntSize
        Offset = Offset + IntSize + k
      END IF

      IF( WriteData ) THEN
        IF( BinaryOutput ) WRITE( VtuUnit ) k
        ! elementtype is fixed to single nodes (==1)
        DO i = 1, NumberOfNodes
          CALL AscBinIntegerWrite( 1 )
        END DO        
        CALL AscBinIntegerWrite( 0, .TRUE. )     
      END IF

      IF( AsciiOutput ) THEN
        WRITE( OutStr,'(A)') lf//'        </DataArray>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF
      IF( WriteXml ) THEN
        WRITE( OutStr,'(A)') '      </Cells>'//lf
        CALL AscBinStrWrite( OutStr ) 
        WRITE( OutStr,'(A)') '    </Piece>'//lf
        CALL AscBinStrWrite( OutStr ) 
        WRITE( OutStr,'(A)') '  </UnstructuredGrid>'//lf
        CALL AscBinStrWrite( OutStr ) 
      END IF


      IF( BinaryOutput ) THEN
        IF( WriteXML ) THEN
          WRITE( OutStr,'(A)') '<AppendedData encoding="raw">'//lf                    
          CALL AscBinStrWrite( OutStr )           
          WRITE( VtuUnit ) '_'
          
          WriteXML = .FALSE.
          WriteData = .TRUE.
          GOTO 100
        ELSE
          WRITE( OutStr,'(A)') lf//'</AppendedData>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF
      
      WRITE( OutStr,'(A)') '</VTKFile>'//lf
      CALL AscBinStrWrite( OutStr ) 
        
      WRITE( OutStr,'(A)') ' '
      CALL AscBinStrWrite( OutStr ) 
      
      CLOSE( VtuUnit )
      
      DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z )

      CALL AscBinWriteFree()
      
    END SUBROUTINE WriteVtuFile
    


    SUBROUTINE WritePvtuFile( VtuFile )
      CHARACTER(LEN=*), INTENT(IN) :: VtuFile
      INTEGER, PARAMETER :: VtuUnit = 58
      TYPE(Variable_t), POINTER :: Var, Solution
      INTEGER :: i,j,k,dofs,Rank,cumn,n,vari,sdofs
      CHARACTER(:), ALLOCATABLE :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
          FieldName2
      LOGICAL :: ScalarsExist, VectorsExist, Found, ComponentVector, ThisOnly
      REAL(KIND=dp), POINTER :: ScalarValues(:), VectorValues(:,:)
      
      
      OPEN( UNIT=VtuUnit, FILE=VtuFile, STATUS='UNKNOWN' )
      
      IF ( LittleEndian() ) THEN
        WRITE( VtuUnit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
      ELSE
        WRITE( VtuUnit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="BigEndian">'
      END IF
      WRITE( VtuUnit,'(A)') '  <PUnstructuredGrid>'
      
      ! nodewise information
      !-------------------------------------
      ScalarFieldName = ListGetString( Params,'Scalar Field 1',ScalarsExist)
      VectorFieldName = ListGetString( Params,'Vector Field 1',VectorsExist)
      IF( ScalarsExist .AND. VectorsExist) THEN
        WRITE( VtuUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)&
            //'" Vectors="'//TRIM(VectorFieldName)//'">'
      ELSE IF( ScalarsExist ) THEN
        WRITE( VtuUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)//'">'
      ELSE IF( VectorsExist ) THEN
        WRITE( VtuUnit,'(A)') '    <PPointData Vectors="'//TRIM(VectorFieldName)//'">'
      END IF


      !-------------------------------------------------------      
      ! Do the scalars
      !-------------------------------------------------------          
      DO Vari = 1, 99
        Txt = 'Scalar Field '//I2S(Vari)
        FieldName = ListGetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        IF( ASSOCIATED( Particles % uvw ) ) THEN
          Solution => VariableGet( Mesh % Variables, TRIM(FieldName), ThisOnly )
          
          IF( .NOT. ASSOCIATED(Solution ) ) THEN
            CALL Warn('WritePvtuFile','Solution not associated!')
          END IF
        ELSE
          IF( FieldName /= 'particle distance' .OR. &
              FieldName /= 'particle time'     .OR. &
              FieldName /= 'particle dt' ) THEN
            CALL Warn('WriteVtuXMLFile','Nonexistent variable: '//TRIM(FieldName))
            CYCLE
          END IF
        END IF
        
        IF( AsciiOutput ) THEN
          WRITE( VtuUnit,'(A)') '      <PDataArray type="Float'//I2S(PrecBits)//&
              '" Name="'//TRIM(FieldName)//'" NumberOfComponents="1" format="ascii"/>'    
        ELSE
          WRITE( VtuUnit,'(A)') '      <PDataArray type="Float'//I2S(PrecBits)//&
              '" Name="'//TRIM(FieldName)//'" NumberOfComponents="1" format="appended"/>'    
        END IF

      END DO


      !-------------------------------------------------------      
      ! Do the vectors
      !-------------------------------------------------------          

      DO Vari = 1, 99
        Txt = 'Vector Field '//I2S(Vari)
        FieldName = ListGetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        IF( ASSOCIATED( Particles % uvw ) ) THEN
          Solution => VariableGet( Mesh % Variables, TRIM(FieldName), ThisOnly )
          ComponentVector = .FALSE.
          
          IF(.NOT. ASSOCIATED(Solution)) THEN
            Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly )
            IF( ASSOCIATED(Solution)) THEN 
              ComponentVector = .TRUE.
            ELSE
              CALL Warn('WriteVtuXMLFile','Nonexistent variable: '//TRIM(FieldName))
              CYCLE
            END IF
          END IF
          
          IF( ASSOCIATED(Solution % EigenVectors)) THEN
            CALL Warn('WritePvtuFile','Do the eigen values')
          END IF
          
          dofs = Solution % DOFs
          IF( ComponentVector ) THEN
            Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly )
            IF( ASSOCIATED(Solution)) dofs = 2
            Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly )
            IF( ASSOCIATED(Solution)) dofs = 3
          END IF
        ELSE
          dofs = 3
          IF( FieldName /= 'velocity' .AND. FieldName /= 'force') THEN
            CALL Warn('WriteVtuXMLFile','Nonexistent variable: '//TRIM(FieldName))
            CYCLE
          END IF
        END IF

        sdofs = dofs
        IF( AsciiOutput ) THEN
          WRITE( VtuUnit,'(A,I1,A)') '      <PDataArray type="Float'//I2S(PrecBits)//'" Name="&
              '//TRIM(FieldName)//'" NumberOfComponents="',sdofs,'" format="ascii"/>'    
        ELSE
          WRITE( VtuUnit,'(A,I1,A)') '      <PDataArray type="Float'//I2S(PrecBits)//'" Name="&
              '//TRIM(FieldName)//'" NumberOfComponents="',sdofs,'" format="appended"/>'    
        END  IF
      END DO

      IF ( ScalarsExist .OR. VectorsExist) THEN
        WRITE( VtuUnit,'(A)') '    </PPointData>'
      END IF
    
      ! Coordinates of each point
      !-------------------------------------
      WRITE( VtuUnit,'(A)') '    <PPoints>'
      IF( AsciiOutput ) THEN
        WRITE( VtuUnit,'(A)') '      <DataArray type="Float'//I2S(PrecBits)//&
            '" NumberOfComponents="3" format="ascii"/>'    
      ELSE
        WRITE( VtuUnit,'(A)') '      <DataArray type="Float'//I2S(PrecBits)//&
        '" NumberOfComponents="3" format="appended"/>'    
      END IF
      WRITE( VtuUnit,'(A)') '    </PPoints>' 
      
      DO i=1,Partitions
        IF( iTime < 10000 ) THEN
          WRITE( VtuUnit,'(A,I4.4,A,I4.4,A)' ) '    <Piece Source="'//&
              TRIM(FilePrefix),i,"par",iTime,'.vtu"/>'
        ELSE
          WRITE( VtuUnit,'(A,I4.4,A,I0,A)' ) '    <Piece Source="'//&
              TRIM(FilePrefix),i,"par",iTime,'.vtu"/>'
        END IF
      END DO
      
      WRITE( VtuUnit,'(A)') '  </PUnstructuredGrid>'
      WRITE( VtuUnit,'(A)') '</VTKFile>'
      
      CLOSE( VtuUnit )
      
    END SUBROUTINE WritePvtuFile
   

!------------------------------------------------------------------------------
  END SUBROUTINE ParticleOutputVtu
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Writes data out in XML VTK ImageData format (VTI) which assumes a uniform grid where
!> the position of each point is defined by the origin and the grid density. 
!> Also binary output and single precision therein is supported. 
!------------------------------------------------------------------------------
  SUBROUTINE ParticleOutputVti( Particles, GridExtent, GridOrigin, GridDx, GridIndex )
!------------------------------------------------------------------------------
    USE AscBinOutputUtils    

    IMPLICIT NONE
    TYPE(Particle_t), POINTER :: Particles  
    INTEGER :: GridExtent(6)
    REAL(KIND=dp) :: GridOrigin(3), GridDx(3)
    INTEGER, POINTER :: GridIndex(:,:,:)

    TYPE(ValueList_t),POINTER :: Params
    INTEGER, SAVE :: nTime = 0
    LOGICAL :: GotIt, Parallel, FixedMeshend
    
    CHARACTER(MAX_NAME_LEN), SAVE :: FilePrefix
    CHARACTER(MAX_NAME_LEN) :: VtiFile, PvtiFile 
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i, j, k, Partitions, Part, ExtCount, FileindexOffSet, iTime
    CHARACTER(:), ALLOCATABLE :: Dir
    REAL(KIND=dp) :: SaveNodeFraction
    LOGICAL :: Found,BinaryOutput,AsciiOutput,SinglePrec,NoFileIndex, &
        Visited = .FALSE.
  
    INTEGER :: NumberOfNodes, ParallelNodes, Dim
    
    
    Params => ListGetSolverParams()
    Mesh => GetMesh()
    
    ExtCount = ListGetInteger( Params,'Output Count',GotIt)
    IF( GotIt ) THEN
      nTime = ExtCount
    ELSE
      nTime = nTime + 1
    END IF
    FileIndexOffset = ListGetInteger( Params,'Fileindex offset',GotIt)
    iTime = nTime + FileIndexOffset

    BinaryOutput = GetLogical( Params,'Binary Output',GotIt)
    AsciiOutput = .NOT. BinaryOutput
    SinglePrec = GetLogical( Params,'Single Precision',GotIt) 
    NoFileindex = GetLogical( Params,'No Fileindex',GotIt)

    IF ( nTime == 1 ) THEN
      FilePrefix = ListGetString( Params,'Filename Prefix')
      CALL Info('ParticleOutputVti','Saving in ImageData VTK XML format to file: ' &
	//TRIM(FilePrefix)//'.vti')
    END IF
    
    Partitions = ParEnv % PEs
    Part = ParEnv % MyPE
    Parallel = (Partitions > 1) .OR. ListGetLogical(Params,'Enforce Parallel format',GotIt)
    
    Dim = Particles % dim
    
    NumberOfNodes = Particles % NumberOfParticles
    
    
    IF (LEN_TRIM(Mesh % Name) > 0 ) THEN
      Dir = TRIM(Mesh % Name) // "/"
    ELSE
      Dir = "./"
    END IF
    
    IF(Parallel .AND. Part == 0) THEN
      CALL Warn('WriteVtiFile','VTK ImageFile not yet in parallel')
      IF(.FALSE.) THEN
        WRITE( PvtiFile,'(A,A,I4.4,".pvti")' ) TRIM(Dir),TRIM(FilePrefix),iTime
        CALL WritePvtiFile( PvtiFile )
      END IF
    END IF
    
    IF ( Parallel ) THEN
      IF( NoFileindex ) THEN
        WRITE( VtiFile,'(A,A,I4.4,A,".vti")' ) TRIM(Dir),TRIM(FilePrefix),Part+1,"par"
      ELSE
        WRITE( VtiFile,'(A,A,I4.4,A,I4.4,".vti")' ) TRIM(Dir),TRIM(FilePrefix),Part+1,"par",&
            iTime
      END IF
    ELSE
      IF( NoFileIndex ) THEN
        WRITE( VtiFile,'(A,A,".vti")' ) TRIM(Dir),TRIM(FilePrefix)
      ELSE
        WRITE( VtiFile,'(A,A,I4.4,".vti")' ) TRIM(Dir),TRIM(FilePrefix),iTime
      END IF
    END IF

    CALL WriteVtiFile( VtiFile )


  CONTAINS


    SUBROUTINE WriteVtiFile( VtiFile )
      CHARACTER(LEN=*), INTENT(IN) :: VtiFile
      INTEGER, PARAMETER :: VtiUnit = 58
      TYPE(Variable_t), POINTER :: Var, Solution, Solution2
      INTEGER :: i,j,k,l,dofs,Rank,cumn,n,vari,sdofs,ind,IsVector,IsAppend,GridPoints,Offset
      CHARACTER(:), ALLOCATABLE :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
          FieldName2, BaseString
      CHARACTER(LEN=1024) :: OutStr
      CHARACTER :: lf
      LOGICAL :: ScalarsExist, VectorsExist, Found, ParticleMode, ComponentVector, &
          ComplementExists, ThisOnly, Stat, WriteData, WriteXML
      INTEGER, POINTER :: Perm(:), Perm2(:), Indexes(:)
      INTEGER, ALLOCATABLE :: ElemInd(:),ElemInd2(:)
      REAL(KIND=dp), POINTER :: ScalarValues(:), VectorValues(:,:),Values(:),Values2(:),&
          Values3(:),Basis(:)
      REAL(KIND=dp) :: x,y,z,u,v,w,DetJ,val
      REAL :: fvalue
      TYPE(Nodes_t) :: Nodes      
      TYPE(Element_t), POINTER :: Element


      ! Initialize the auxiliary module for buffered writing
      !--------------------------------------------------------------
      CALL AscBinWriteInit( AsciiOutput, SinglePrec, VtiUnit, NumberOfNodes )
      
      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

      n = Mesh % MaxElementDOFS
      ALLOCATE( ElemInd(n), ElemInd2(n) )
      ThisOnly = .TRUE.

      ! Linefeed character
      !-----------------------------------
      lf = CHAR(10)


      ! This format works both for ascii and binary output
      !-------------------------------------------------------------------------
      OPEN( UNIT=VtiUnit, FILE=VtiFile, FORM = 'unformatted', ACCESS = 'stream' )


      WRITE( OutStr,'(A)') '<?xml version="1.0"?>'//lf
      CALL AscBinStrWrite( OutStr ) 
            
      IF ( LittleEndian() ) THEN
        WRITE( OutStr,'(A)') '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'//lf
      ELSE
        WRITE( OutStr,'(A)') '<VTKFile type="ImageData" version="0.1" byte_order="BigEndian">'//lf
      END IF
      CALL AscBinStrWrite( OutStr ) 

      WRITE( OutStr,'(A,6I4,A,3ES16.7E3,A,3ES16.7E3,A)') &
          '<ImageData WholeExtent = "',GridExtent,&          
          '"  Origin = "',GridOrigin(1:3), &          
          '"  Spacing = "',GridDx(1:3),'">'//lf
      CALL AscBinStrWrite( OutStr )  

      WRITE( OutStr,'(A,6I4,A)') '  <Piece Extent="',GridExtent,'">'//lf
      CALL AscBinStrWrite( OutStr ) 

      
      !---------------------------------------------------------------------
      ! Header information for nodewise data
      !---------------------------------------------------------------------
      ScalarFieldName = ListGetString( Params,'Scalar Field 1',ScalarsExist)
      VectorFieldName = ListGetString( Params,'Vector Field 1',VectorsExist)
      IF( ScalarsExist .AND. VectorsExist) THEN
        WRITE( OutStr,'(A)') '    <PointData Scalars="'//TRIM(ScalarFieldName)&
            //'" Vectors="'//TRIM(VectorFieldName)//'">'//lf
      ELSE IF( ScalarsExist ) THEN
        WRITE( OutStr,'(A)') '    <PointData Scalars="'//TRIM(ScalarFieldName)//'">'//lf
      ELSE IF( VectorsExist ) THEN
        WRITE( OutStr,'(A)') '    <PointData Vectors="'//TRIM(VectorFieldName)//'">'//lf
      ELSE
        CALL Warn('WriteVtiFile','Are there really no scalars or vectors?')
      END IF
      IF( ScalarsExist .OR. VectorsExist ) THEN
        CALL AscBinStrWrite( OutStr ) 
      END IF
      
      Offset = 0
      WriteXML = .TRUE.
      WriteData = AsciiOutput
      GridPoints = 1
      DO i = 1,3
        GridPoints = GridPoints * ( GridExtent(2*i)-GridExtent(2*i-1)+1 )
      END DO

      
  
100   DO IsVector = 0, 1        

        DO Vari = 1, 99

          !---------------------------------------------------------------------
          ! do the scalars
          !--------------------------------- -----------------------------------
          IF( IsVector == 0 ) THEN
            BaseString = 'Scalar Field'
          ELSE
            BaseString = 'Vector Field'
          END IF
          
          Txt = TRIM(BaseString)//' '//I2S(Vari)          
          FieldName = ListGetString( Params, TRIM(Txt), Found )
          IF(.NOT. Found) EXIT
          
          !---------------------------------------------------------------------
          ! Find the variable with the given name in the normal manner 
          !---------------------------------------------------------------------
          Solution => VariableGet( Mesh % Variables,TRIM(FieldName),ThisOnly )
          IF( ASSOCIATED( Solution ) ) THEN
            IF( IsVector == 1 .AND. Solution % Dofs <= 1 ) NULLIFY( Solution ) 
          END IF

          ComponentVector = .FALSE.
          IF( IsVector == 1 ) THEN
            IF(.NOT. ASSOCIATED(Solution)) THEN
              Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 1', ThisOnly )
              ComponentVector = ASSOCIATED( Solution ) 
            END IF
          END IF
          IF( .NOT. ASSOCIATED( Solution ) ) THEN
            CALL Warn('WriteVtiXMLFile','Nonexistent variable: '//TRIM(FieldName))
            CYCLE
          END IF

          Perm => Solution % Perm
          dofs = Solution % DOFs
          Values => Solution % Values
          
          !---------------------------------------------------------------------
          ! Eigenmodes have not yet been implemented
          !---------------------------------------------------------------------
          IF( ASSOCIATED(Solution % EigenVectors)) THEN
            CALL Warn('WriteVtiXMLFile','Do the eigen values')
          END IF

          !---------------------------------------------------------------------
          ! Some vectors are defined by a set of components (either 2 or 3)
          !---------------------------------------------------------------------
          IF( ComponentVector ) THEN
            Solution2 => VariableGet( Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly )
            IF( ASSOCIATED(Solution2)) THEN
              Values2 => Solution2 % Values
              dofs = 2
            END IF
            Solution2 => VariableGet( Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly )
            IF( ASSOCIATED(Solution2)) THEN
              Values3 => Solution2 % Values
              dofs = 3
            END IF
          END IF
          
          !---------------------------------------------------------------------
          ! There may be special complementary variables such as 
          ! displacement & mesh update 
          !---------------------------------------------------------------------
          ComplementExists = .FALSE.
          Txt = TRIM(BaseString)//I2S(Vari)//' Complement'

          FieldName2 = ListGetString( Params, TRIM(Txt), Found )
          IF( Found ) THEN
            Solution2 => VariableGet( Mesh % Variables, &
                TRIM(FieldName2), ThisOnly )
            IF( ASSOCIATED(Solution2)) THEN 
              Values2 => Solution2 % Values
              Perm2 => Solution2 % Perm 
              ComplementExists = .TRUE.
            ELSE
              CALL Warn('WriteVTIFile','Complement does not exist:'//TRIM(FieldName2))
            END IF
          END IF
          
          
          !---------------------------------------------------------------------
          ! Finally save the field values for scalars and vectors
          !---------------------------------------------------------------------
          j = 0
          
          IF( dofs == 1 ) THEN
            sdofs = 1
          ELSE
            sdofs = MAX(dofs,3)
          END IF
          
          IF( WriteXML ) THEN
            WRITE( OutStr,'(A)') '      <DataArray type="Float64" Name="'//TRIM(FieldName)//'"'
            CALL AscBinStrWrite( OutStr ) 

            WRITE( OutStr,'(A,I1,A)') ' NumberOfComponents="',sdofs,'"'          
            CALL AscBinStrWrite( OutStr ) 

            IF( AsciiOutput ) THEN
              WRITE( OutStr,'(A)') ' format="ascii">'//lf
              CALL AscBinStrWrite( OutStr ) 
            ELSE
              WRITE( OutStr,'(A)') ' format="appended"'
              CALL AscBinStrWrite( OutStr ) 
              WRITE( OutStr,'(A,I8,A)') ' offset="',Offset,'"/>'//lf
              CALL AscBinStrWrite( OutStr ) 
            END IF
          END IF
          
          IF( BinaryOutput ) THEN
            IF( SinglePrec ) THEN
              k = GridPoints * KIND( fvalue ) * sdofs
            ELSE
              k = GridPoints * KIND( val ) * sdofs
            END IF
            Offset = Offset + KIND(i) + k
            IF( WriteData ) WRITE( VtiUnit ) k
          END IF


          IF( WriteData ) THEN
            DO k = 1,GridExtent(6)-GridExtent(5)+1
              DO j = 1,GridExtent(4)-GridExtent(3)+1
                DO i = 1,GridExtent(2)-GridExtent(1)+1
                  
                  ind = GridIndex( i, j, k ) 
                  
                  IF( ind == 0 ) THEN
                    DO l=1,sdofs                  
                      IF( AsciiOutput ) THEN
                        WRITE( OutStr,'(A)') ' 0.0'  
                        CALL AscBinStrWrite( OutStr )
                      ELSE IF( SinglePrec ) THEN
                        fvalue = 0.0
                        WRITE( VtiUnit ) fvalue                          
                      ELSE
                        val = 0.0_dp
                        WRITE( VtiUnit ) val                          
                      END IF
                    END DO
                  ELSE
                    Element => Mesh % Elements( Particles % ElementIndex(ind) )            
                    Indexes => Element % NodeIndexes
                    n = Element % TYPE % NumberOfNodes
                    
                    Nodes % x(1:n) = Mesh % Nodes % x( Indexes ) 
                    Nodes % y(1:n) = Mesh % Nodes % y( Indexes ) 
                    Nodes % z(1:n) = Mesh % Nodes % z( Indexes ) 
                    
                    u = Particles % uvw(ind,1)
                    v = Particles % uvw(ind,2)
                    IF( dim == 3 ) THEN
                      w = Particles % uvw(ind,3)
                    ELSE
                      w = 0.0_dp
                    END IF
                    
                    stat = ElementInfo( Element,Nodes,u,v,w,detJ,Basis)
                    
                    IF( Solution % TYPE == Variable_on_nodes_on_elements ) THEN
                      ElemInd(1:n) = Perm( Element % DGIndexes(1:n) )
                      IF( ComplementExists ) THEN
                        ElemInd2(1:n) = Perm2( Element % DGIndexes(1:n) )
                      END IF
                    ELSE 
                      ElemInd(1:n) = Perm( Indexes(1:n) )
                      IF( ComplementExists ) THEN
                        ElemInd2(1:n) = Perm2( Indexes(1:n) )
                      END IF
                    END IF

                    DO l=1,sdofs                  
                      val = 0.0_dp
                      IF( l <= dofs ) THEN
                        IF( ALL( ElemInd(1:n) > 0 ) ) THEN
                          IF( ComponentVector ) THEN
                            IF( l == 1 ) THEN
                              val = SUM( Basis(1:n) * Values(ElemInd(1:n)) )
                            ELSE IF( l == 2 ) THEN
                              val = SUM( Basis(1:n) * Values2(ElemInd(1:n)) )
                            ELSE
                              val = SUM( Basis(1:n) * Values3(ElemInd(1:n)) )
                            END IF
                          ELSE 
                            val = SUM( Basis(1:n) * Values(dofs*(ElemInd(1:n)-1)+l) )
                          END IF
                        ELSE IF( ComplementExists ) THEN
                          IF( ALL( ElemInd2(1:n) > 0 ) ) THEN
                            val = SUM( Basis(1:n) * Values2(dofs*(ElemInd2(1:n)-1)+l) )
                          END IF
                        END IF
                      END IF
                      
                      IF( AsciiOutput ) THEN
                        IF( ABS( val ) < TINY( val ) ) THEN
                          WRITE( OutStr,'(A)') ' 0.0'  
                        ELSE 
                          WRITE( OutStr,'(ES16.7E3)') val
                        END IF
                        CALL AscBinStrWrite( OutStr )
                      ELSE IF( SinglePrec ) THEN
                        fvalue = val
                        WRITE( VtiUnit ) fvalue
                      ELSE
                        WRITE( VtiUnit ) val
                      END IF

                    END DO
                  END IF
                  
                END DO ! i
              END DO ! j
            END DO ! k
          END IF

          IF( AsciiOutput ) THEN
            WRITE( OutStr,'(A)') lf//'      </DataArray>'//lf
            CALL AscBinStrWrite( OutStr ) 
          END IF

        END DO
      END DO
        

      IF( WriteXML ) THEN
        IF( ScalarsExist .OR. VectorsExist) THEN
          WRITE( OutStr,'(A)') '    </PointData>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
        
        WRITE( OutStr,'(A)') '    <CellData>'//lf
        CALL AscBinStrWrite( OutStr ) 

        WRITE( OutStr,'(A)') '    </CellData>'//lf
        CALL AscBinStrWrite( OutStr ) 
       
        WRITE( OutStr,'(A)') '  </Piece>'//lf
        CALL AscBinStrWrite( OutStr ) 

        WRITE( OutStr,'(A)') '</ImageData>'//lf
        CALL AscBinStrWrite( OutStr )         
      END IF

      IF( BinaryOutput ) THEN
        IF( WriteXML ) THEN
          WRITE( OutStr,'(A)') '<AppendedData encoding="raw">'//lf                    
          CALL AscBinStrWrite( OutStr )           
          WRITE( VtiUnit ) '_'
          
          WriteXML = .FALSE.
          WriteData = .TRUE.
          GOTO 100
        ELSE
          WRITE( OutStr,'(A)') lf//'</AppendedData>'//lf
          CALL AscBinStrWrite( OutStr ) 
        END IF
      END IF


      WRITE( OutStr,'(A)') '</VTKFile>'//lf
      CALL AscBinStrWrite( OutStr ) 
      
      CLOSE( VtiUnit )
      
      DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z, ElemInd, ElemInd2 )

      CALL AscBinWriteFree()

      
    END SUBROUTINE WriteVtiFile
      
    
    
    SUBROUTINE WritePvtiFile( VtiFile )
      CHARACTER(LEN=*), INTENT(IN) :: VtiFile
      INTEGER, PARAMETER :: VtiUnit = 58
      TYPE(Variable_t), POINTER :: Var, Solution
      INTEGER :: i,j,k,dofs,Rank,cumn,n,vari,sdofs
      CHARACTER(:), ALLOCATABLE :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
          FieldName2
      LOGICAL :: ScalarsExist, VectorsExist, Found, ComponentVector, ThisOnly
      REAL(KIND=dp), POINTER :: ScalarValues(:), VectorValues(:,:)
      
      
      OPEN( UNIT=VtiUnit, FILE=VtiFile, STATUS='UNKNOWN' )
      
      IF ( LittleEndian() ) THEN
        WRITE( VtiUnit,'(A)') '<VTKFile type="PImageData" version="0.1" byte_order="LittleEndian">'
      ELSE
        WRITE( VtiUnit,'(A)') '<VTKFile type="PImageData" version="0.1" byte_order="BigEndian">'
      END IF
      WRITE( VtiUnit,'(A)') '  <PImageData>'
      
      ! nodewise information
      !-------------------------------------
      ScalarFieldName = ListGetString( Params,'Scalar Field 1',ScalarsExist)
      VectorFieldName = ListGetString( Params,'Vector Field 1',VectorsExist)
      IF( ScalarsExist .AND. VectorsExist) THEN
        WRITE( VtiUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)&
            //'" Vectors="'//TRIM(VectorFieldName)//'">'
      ELSE IF( ScalarsExist ) THEN
        WRITE( VtiUnit,'(A)') '    <PPointData Scalars="'//TRIM(ScalarFieldName)//'">'
      ELSE IF( VectorsExist ) THEN
        WRITE( VtiUnit,'(A)') '    <PPointData Vectors="'//TRIM(VectorFieldName)//'">'
      END IF
      
      
      !-------------------------------------------------------      
      ! Do the scalars
      !-------------------------------------------------------          
      DO Vari = 1, 99
        Txt = 'Scalar Field '//I2S(Vari)
        FieldName = ListGetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        Solution => VariableGet( Mesh % Variables, TRIM(FieldName), ThisOnly )
        
        IF( .NOT. ASSOCIATED(Solution ) ) THEN
          CALL Warn('WritePvtiFile','Solution not associated!')
        END IF
        
        WRITE( VtiUnit,'(A)') '      <PDataArray type="Float64" Name="'//TRIM(FieldName)&
            //'" NumberOfComponents="1" format="ascii"/>'    
      END DO
      
      
      !-------------------------------------------------------      
      ! Do the vectors
      !-------------------------------------------------------          
      
      DO Vari = 1, 99
        Txt = 'Vector Field '//I2S(Vari)
        FieldName = ListGetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        Solution => VariableGet( Mesh % Variables, TRIM(FieldName), ThisOnly )
        ComponentVector = .FALSE.
        
        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly )
          IF( ASSOCIATED(Solution)) THEN 
            ComponentVector = .TRUE.
          ELSE
            CALL Warn('WriteVtiXMLFile', 'Nonexistent variable: '//TRIM(FieldName))
            CYCLE
          END IF
        END IF
        
        IF( ASSOCIATED(Solution % EigenVectors)) THEN
          CALL Warn('WritePvtiFile','Do the eigen values')
        END IF
        
        dofs = Solution % DOFs
        IF( ComponentVector ) THEN
          Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 2',ThisOnly )
          IF( ASSOCIATED(Solution)) dofs = 2
          Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 3',ThisOnly )
          IF( ASSOCIATED(Solution)) dofs = 3
        END IF
        
        sdofs = dofs
        WRITE( VtiUnit,'(A,I1,A)') '      <PDataArray type="Float64" Name="'//TRIM(FieldName)&
            //'" NumberOfComponents="',sdofs,'" format="ascii"/>'    
      END DO
      
      WRITE( VtiUnit,'(A)') '    </PPointData>'
      
      ! Coordinates of each point
      !-------------------------------------
      WRITE( VtiUnit,'(A)') '    <PPoints>'
      WRITE( VtiUnit,'(A)') '      <DataArray type="Float64" NumberOfComponents="3" format="ascii"/>'    
      WRITE( VtiUnit,'(A)') '    </PPoints>' 
      
      DO i=1,Partitions
        IF( NoFileindex ) THEN
          WRITE( VtiUnit,'(A,I4.4,A,A)' ) '    <Piece Source="'//&
              TRIM(FilePrefix),i,"par",'.vti"/>'
        ELSE
          WRITE( VtiUnit,'(A,I4.4,A,I4.4,A)' ) '    <Piece Source="'//&
              TRIM(FilePrefix),i,"par",iTime,'.vti"/>'
        END IF
      END DO
      
      WRITE( VtiUnit,'(A)') '  </PImageData>'
      WRITE( VtiUnit,'(A)') '</VTKFile>'
      
      CLOSE( VtiUnit )
      
    END SUBROUTINE WritePvtiFile
   

!------------------------------------------------------------------------------
  END SUBROUTINE ParticleOutputVti
!------------------------------------------------------------------------------


 !------------------------------------------------------------------------
 !> Write particles to an external file in various formats.
 !-------------------------------------------------------------------------

  SUBROUTINE SaveParticleData( Model,Solver,dt,TransientSimulation )
    
    
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    
    TYPE(Particle_t), POINTER :: Particles
    LOGICAL :: Visited = .FALSE.
    TYPE(ValueList_t), POINTER :: Params     
    CHARACTER(LEN=MAX_NAME_LEN) :: FileFormat
    LOGICAL :: VtuFormat, TableFormat, GmshFormat, AnyFormat, Found
    
    SAVE :: TableFormat, VtuFormat, Visited

    
    Particles => GlobalParticles
    Params => ListGetSolverParams()
    
    TableFormat = ListGetLogical( Params,'Table Format',Found)
    GmshFormat = ListGetLogical( Params,'Gmsh Format',Found)
    VtuFormat = ListGetLogical( Params,'Vtu Format',Found)
    FileFormat = ListGetString( Params,'Output Format',Found) 
    IF( Found ) THEN
      IF( FileFormat == 'gmsh') GmshFormat = .TRUE.
      IF( FileFormat == 'vtu') VtuFormat = .TRUE.
      IF( FileFormat == 'table') TableFormat = .TRUE.
    END IF
    
    AnyFormat = TableFormat .OR. VtuFormat .OR. GmshFormat
    IF( .NOT. AnyFormat ) THEN
      CALL Warn('SaveParticleData','No active file format given!')
      RETURN
    END IF
    
    IF(.NOT. ListCheckPresent( Params,'Filename Prefix') ) THEN
      CALL ListAddString( Params,'Filename Prefix','particles')
    END IF
    
    IF( TableFormat ) CALL ParticleOutputTable( Particles )
    IF( GmshFormat ) CALL ParticleOutputGmsh( Particles ) 
    IF( VtuFormat ) CALL ParticleOutputVtu( Particles ) 


  END SUBROUTINE SaveParticleData
  

!------------------------------------------------------------------------------
END MODULE ParticleUtils
!------------------------------------------------------------------------------

!> \}

