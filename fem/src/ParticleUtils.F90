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
    k = NINT( ParallelReduction( 1.0_dp * NoParticles ) )
    WRITE(Message,'(A,T18,I0)') 'Total: ',k
    CALL Info('ParticleStatusCount',Message,Level=8)
    DO i=1,PARTICLE_GHOST
      j = StatusCount(i)
      k = NINT( ParallelReduction( 1.0_dp * j ) )
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
    
    i = NINT( ParallelReduction( 1.0_dp * i ) )
    j = NINT( ParallelReduction( 1.0_dp * j ) )

    CALL Info('MarkInternalElements','Internal Elements: '//TRIM(I2S(i)),Level=8 )
    CALL Info('MarkInternalElements','Interface Elements: '//TRIM(I2S(j)),Level=8 )    
    
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
        //TRIM(I2S(n2)),Level=12)

    IF(n1 == 0 ) THEN
      CALL Info('DeleteLostParticles','No particles need to be deleted',Level=12)
      RETURN
    ELSE
      CALL Info('DeleteLostParticles','First particle with changed permutation: '&
          //TRIM(I2S(n1)),Level=12)      
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
        //TRIM(I2S(ReleaseCount)),Level=10)
    
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
    
    INTEGER, POINTER :: Neighbours(:)
    LOGICAL, POINTER :: FaceInterface(:), IsNeighbour(:)

    INTEGER :: q
    LOGICAL, ALLOCATABLE :: Failed(:)
    
    TYPE ExchgInfo_t
      INTEGER :: n=0
      INTEGER, ALLOCATABLE :: Gindex(:), Particles(:)
    END TYPE ExchgInfo_t
    
    REAL(KIND=dp), ALLOCATABLE :: Buf(:)
    INTEGER, ALLOCATABLE :: BufInt(:)
    TYPE(ExchgInfo_t), ALLOCATABLE :: ExcInfo(:)
    !---------------------------------------------------------
    
    nReceived = 0
    IF( ParEnv% PEs == 1 ) RETURN
    
    CALL Info('ChangeParticlePartition','Sending particles among partitions',Level=10)

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

    CALL Info('ChangeParticlePartition','Number of active partitions: '&
        //TRIM(I2S(NoPartitions)),Level=12)
    
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
        IF ( Face % BoundaryInfo % Constraint > 0 ) &
            CALL Warn("ChangeParticlePartition", "is this a BC after all?")
        
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
   
    n = SUM( ExcInfo(1:NoPartitions) % n )
    CALL Info('ChangeParticlePartition','Number of particles to send: '&
        //TRIM(I2S(n)),Level=10)
    
    CALL MPI_ALLREDUCE( n, nSent, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )
    IF ( nSent == 0 ) THEN
      CALL Info('ChangeParticlePartition','No particles needs to be sent',Level=10)
      DEALLOCATE(ExcInfo, Perm, Neigh)
      RETURN
    ELSE
      CALL Info('ChangeParticlePartition','Global number of particles to sent: '&
          //TRIM(I2S(nSent)),Level=10)      
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
    CALL Info('ChangeParticlePartition','Number of particles to recieve: '&
        //TRIM(I2S(n)),Level=10)
   
    CALL MPI_ALLREDUCE( n, nReceived, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )
    IF ( nReceived==0 ) THEN
      CALL Info('ChangeParticlePartition','No particles needs to be received',Level=10)
      DEALLOCATE(Recv_Parts, Requests, ExcInfo, Perm, Neigh)
      RETURN
    ELSE
      CALL Info('ChangeParticlePartition','Global number of particles to recieve: '&
          //TRIM(I2S(nReceived)),Level=10)      
    END IF
   
    n = SUM( ExcInfo(1:NoPartitions) % n )
    CALL Info('ChangeParticlePartition','Total number of particles to sent: '&
        //TRIM(I2S(n)),Level=10)


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
          CALL Warn( 'ChangeParticlePartition', 'Neighbouring partition not found?')
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

    CALL Info('ChangeParticlePartition','Collected particles from partitions: '&
        //TRIM(I2S(n)),Level=12)

    ncomp = dim ! coordinate
    IF( ASSOCIATED( Particles % Velocity ) ) ncomp = ncomp + dim
    IF( ASSOCIATED( Particles % PrevCoordinate ) ) ncomp = ncomp + dim
    IF( ASSOCIATED( Particles % PrevVelocity ) ) ncomp = ncomp + dim
    IF( ASSOCIATED( Particles % Force ) ) ncomp = ncomp + dim
    Var => Particles % Variables
    DO WHILE( ASSOCIATED(Var) )
      ncomp = ncomp + Var % Dofs 
      IF( Var % Dofs /= 1 ) CALL Warn('ChangeParticlePartition','Implement for vectors!')
      Var => Var % Next 
    END DO

    ncompInt = 0
    IF ( ASSOCIATED(Particles % NodeIndex) ) ncompInt = ncompInt + 1
    IF ( ASSOCIATED(Particles % Partition) ) ncompInt = ncompInt + 1
    ! status, elementindex & closestnode are recomputed, and hence not communicated

    CALL Info('ChangeParticlePartition','Transferring real entries between particles: '&
        //TRIM(I2S(ncomp)),Level=12)
    CALL Info('ChangeParticlePartition','Transferring integer entries between particles: '&
        //TRIM(I2S(ncompInt)),Level=12)


    n = 2*(n + 2*ncomp + MPI_BSEND_OVERHEAD*2*NoPartitions)
    CALL CheckBuffer(n)
    
    CALL Info('ChangeParticlePartition','Size of data buffer: ' &
        //TRIM(I2S(n)),Level=12)

    ! Send particles:
    ! ---------------
    CALL Info('ChangeParticlePartition','Now sending particle data',Level=14)
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
    CALL Info('ChangeParticlePartition','Now recieving particle data',Level=14)
    
    n = SUM(Recv_Parts)

    CALL DeleteLostParticles(Particles)
    IF ( Particles % NumberOfParticles+n > Particles % MaxNumberOfParticles ) THEN
      CALL IncreaseParticles( Particles, Particles % NumberOfParticles + 2*n - &
                    Particles % MaxNumberOfParticles )
    ELSE
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
      
      Particles % NumberOfParticles = Particles % NumberOfParticles + COUNT(.NOT.Failed)
      DEALLOCATE(Indexes, Failed)
    END DO
    
    DEALLOCATE(Recv_Parts, Neigh, Requests)
    CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )

    CALL Info('ChangeParticlePartition','Information exchange done',Level=10)

    
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
  

  !---------------------------------------------------------
  !> Subroutine for advecting back to given partition and 
  !> node index.
  !---------------------------------------------------------
  SUBROUTINE ParticleAdvectParallel(Particles, SentField, RecvField, dofs ) 
    !---------------------------------------------------------
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp), POINTER :: SentField(:), RecvField(:)
    INTEGER :: Dofs
    !---------------------------------------------------------
    INTEGER i,j,k,l,m,n,dim,NoPartitions,NoParticles,part, &
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
    IF( Dofs /= 1 ) THEN
      CALL Fatal('ParticleAdvectParallel','Implement for more than one dof!')
    END IF

    ! First take the components of the field that lie in the present partition. 
    !-------------------------------------------------------------------------
    DO i=1,Particles % NumberOfParticles
      Part = Particles % Partition(i)
      IF( Part-1 == ParEnv % MyPe ) THEN
        j = Particles % NodeIndex(i)
        IF ( j==0 ) CYCLE
        RecvField(j) = SentField(i)
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
      CALL Info('ParticleAdvectParallel','Invalid partition in particles: '//TRIM(I2S(nerr)))
    END IF
    
    n = SUM( SentParts )
    CALL Info('ParticleAdvectParallel','Local particles to be sent: '//TRIM(I2S(n)),Level=12)
    
    CALL MPI_ALLREDUCE( n, nSent, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )

    IF( nSent > 0 ) THEN
      CALL Info('ParticleAdvectParallel','Global particles to be sent: '&
          //TRIM(I2S(nSent)),Level=12)
    ELSE
      ! If nobody is sending any particles then there can be no need to reveive particles either
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
    CALL Info('ParticleAdvectParallel','Particles to be recieved: '//TRIM(I2S(n)),Level=12)

    CALL MPI_ALLREDUCE( n, nReceived, 1, MPI_INTEGER, &
        MPI_SUM, ELMER_COMM_WORLD, ierr )

    IF ( nReceived==0 ) THEN
      DEALLOCATE(RecvParts, SentParts, Requests )
      RETURN
    END IF

    CALL Info('ParticleAdvectParallel','Total number of particles to be recieved: '&
        //TRIM(I2S(nReceived)),Level=12)    

    n = 2*(n + MPI_BSEND_OVERHEAD*2*NoPartitions)
    CALL CheckBuffer(n)

    CALL Info('ParticleAdvectParallel','Buffer size for sending and recieving: ' &
        //TRIM(I2S(n)),Level=14)    


    ! Allocate sent and recieve buffers based on the maximum needed size.
    !--------------------------------------------------------------------
    n = MAXVAL( SentParts )
    ALLOCATE( SentReal(n), SentInt(n) )
    CALL Info('ParticleAdvectParallel','Allocating sent buffer of size: '&
        //TRIM(I2S(n)),Level=18)

    n = MAXVAL( RecvParts ) 
    ALLOCATE( RecvReal(n), RecvInt(n) )
    CALL Info('ParticleAdvectParallel','Allocating recieve buffer of size: '&
        //TRIM(I2S(n)),Level=18)
   
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
          SentReal(m) = SentField(l)
          SentInt(m) = Particles % NodeIndex(l)
        END IF
      END DO

      CALL MPI_BSEND( SentInt, m, MPI_INTEGER, j-1, &
          1001, ELMER_COMM_WORLD, ierr )

      CALL MPI_BSEND( SentReal, m, MPI_DOUBLE_PRECISION, j-1, &
          1002, ELMER_COMM_WORLD, ierr )
    END DO


    ! Recv particles:
    ! ---------------   
    CALL Info('ParticleAdvectParallel','Now recieving field',Level=14)

    nerr = 0 
    DO j=1,NoPartitions

      IF( j-1 == ParEnv % MyPe ) CYCLE

      m = RecvParts(j)
      IF ( m == 0 ) CYCLE
      
      CALL MPI_RECV( RecvInt, m, MPI_INTEGER, j-1, &
          1001, ELMER_COMM_WORLD, status, ierr )

      CALL MPI_RECV( RecvReal, m, MPI_DOUBLE_PRECISION, j-1, &
          1002, ELMER_COMM_WORLD, status, ierr )
      
      DO l=1,m
        k = RecvInt(l)
        IF( k <=0 .OR. k > SIZE( RecvField ) ) THEN
          nerr = nerr + 1
          CYCLE
        END IF
        RecvField(k) = RecvReal(l)
      END DO
    END DO
        
    IF( nerr > 0 ) THEN
      CALL Info('ParticleAdvectParallel','Invalid recieved index in particles: '//TRIM(I2S(nerr)))
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
    REAL(KIND=dp), POINTER :: TargetVector(:,:)
    INTEGER, POINTER :: Status(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: DataName

    
    MeanCoord = 0.0_dp
    AbsCoord = 0.0_dp
    VarCoord = 0.0_dp
    MinCoord = HUGE( MinCoord )
    MaxCoord = -HUGE( MaxCoord )
    
    Cnt = 0
    NoParticles =  Particles % NumberOfParticles
    dim = Particles % dim
    Coord = 0.0_dp
    
    IF( DerOrder == -1 ) THEN
      TargetVector => Particles % PrevCoordinate
      DataName = 'previous coordinate values'
    ELSE IF( DerOrder == 0 ) THEN
      TargetVector => Particles % Coordinate
      DataName = 'current coordinate values'      
    ELSE IF( DerOrder == 1 ) THEN
      TargetVector => Particles % Velocity
      DataName = 'current velocity values'
    ELSE IF( DerOrder == 2 ) THEN
      TargetVector => Particles % Force
      DataName = 'current force values'
    ELSE
      CALL Fatal('ParticleStatistics','Unknown value for DerOrder!')
    END IF
    
    Status => Particles % Status
    
    DO i=1,NoParticles
      IF( Status(i) >= PARTICLE_LOST ) CYCLE
      IF( Status(i) < PARTICLE_INITIATED ) CYCLE
      
      Coord(1:dim) = TargetVector(i,1:dim)
      
      MeanCoord = MeanCoord + Coord
      AbsCoord = AbsCoord + ABS( Coord )
      VarCoord = VarCoord + Coord**2
      DO j=1,dim
        MinCoord(j) = MIN( MinCoord(j), Coord(j) )
        MaxCoord(j) = MAX( MaxCoord(j), Coord(j) )
      END DO
      Cnt = Cnt + 1
    END DO
    
    TotParticles = NINT( ParallelReduction( 1.0_dp * Cnt ) )
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
      TotNoParticles =  NINT( ParallelReduction( 1.0_dp * Particles % NumberOfParticles ) )
      TotParticleStepsTaken = NINT( ParallelReduction( 1.0_dp * ParticleStepsTaken) )
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
      ParallelParticles = NINT( ParallelReduction( 1.0_dp * Cnt ) )
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
      n = Element % TYPE % NumberOfNodes

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
    NoElems = NINT( ParallelReduction( 1.0_dp * NoElems ) )
    
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
          
    Visited = .TRUE.
    
  END FUNCTION CharacteristicElementSize


  
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
  SUBROUTINE InitializeParticles( Particles, InitParticles, AppendParticles, Group, SaveOrigin ) 
    
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: InitParticles
    LOGICAL, OPTIONAL :: AppendParticles
    INTEGER, OPTIONAL :: Group
    LOGICAL, OPTIONAL :: SaveOrigin
    
    TYPE(ValueList_t), POINTER :: Params, BodyForce 
    TYPE(Variable_t), POINTER :: Var, AdvVar
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t) :: Nodes
    INTEGER :: Offset, NewParticles,LastParticle,NoElements
    INTEGER :: dim, ElementIndex, body_id, bf_id
    REAL(KIND=dp), POINTER :: rWork(:,:),Coordinate(:,:), Velocity(:,:)
    REAL(KIND=dp) :: Velo(3), Coord(3), Center(3), CenterVelo(3), time0, dist
    CHARACTER(LEN=MAX_NAME_LEN) :: InitMethod
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
    LOGICAL :: CheckForSize, Parallel, SaveParticleOrigin
    LOGICAL, POINTER :: DoneParticle(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, str
    
    SAVE Nodes


    Mesh => GetMesh()
    Params => ListGetSolverParams()
    dim = Particles % Dim
    Parallel = ( ParEnv % PEs > 1 )

    IF( PRESENT( SaveOrigin ) ) THEN
      SaveParticleOrigin = SaveOrigin
    ELSE
      SaveParticleOrigin = .FALSE.
    END IF
    
    !------------------------------------------------------------------------
    ! Initialize the timestepping strategy stuff
    !-------------------------------------------------------------------------
    
    InitMethod = ListGetString( Params,'Coordinate Initialization Method',gotIt ) 
    Particles % RK2 = ListGetLogical( Params,'Runge Kutta', GotIt )
    IF( Particles % RK2 .AND. ParEnv % PEs > 1 ) THEN
      CALL Warn('InitializeParticles','> Runge Kutta < integration might not work in parallel')
    END IF
    
    Particles % DtConstant = ListGetLogical( Params,'Particle Dt Constant',GotIt )
    IF(.NOT. GotIt) Particles % DtConstant = .TRUE.

    IF( ListGetLogical( Params,'Particle Dt Negative',GotIt ) ) THEN
      Particles % DtSign = -1
    END IF

    
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
      Var => VariableGet( Mesh % Variables, TRIM(VariableName) )
      IF( .NOT. ASSOCIATED( Var ) ) THEN
        CALL Fatal('InitializeParticles','Mask / Condition variable does not exist!')
      END IF

      MaskPerm => Var % Perm
      MaskVal => Var % Values

      IF(.NOT. ( ASSOCIATED( MaskPerm ) .AND. ASSOCIATED(MaskVal)) ) THEN
        CALL Warn('InitializeParticles','Initialization variable does not exist?')
      ELSE IF( MAXVAL( MaskPerm ) == 0 ) THEN
        CALL Warn('InitializeParticles','Initialization variable of size zero?')
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

          PRINT *,'Total nodes vs. masked',Mesh % NumberOfNodes,nonodes
        ELSE IF( InitMethod == 'elemental') THEN
          ALLOCATE( InvPerm( MAX( Mesh % NumberOfBulkElements, Mesh % NumberOfBoundaryElements ) ) ) 
          InvPerm = 0

          j = 0
          DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
            CurrentElement => Mesh % Elements(i)
            NodeIndexes =>  CurrentElement % NodeIndexes
            n = CurrentElement % TYPE % NumberOfNodes

            IF( i == Mesh % NumberOfBulkElements ) THEN
              IF( j > 0 ) EXIT
            END IF

            IF( ANY( MaskPerm( NodeIndexes ) == 0 ) ) CYCLE

            IF( RequirePositivity ) THEN
              meanval = SUM( MaskVal( MaskPerm( NodeIndexes ) ) ) 
              IF( meanval < 0.0_dp ) CYCLE
            END IF

            ! If some of bulk elements have been found avtive
            j = j + 1
            InvPerm(j) = i

          END DO
          noelements = j

          PRINT *,'Total elements vs. masked',Mesh % NumberOfBulkElements,noelements
        END IF
      END IF
    ELSE
      nonodes = Mesh % NumberOfNodes
      noelements = Mesh % NumberOfBulkElements
    END IF
    
    GotWeight = ListCheckPresentAnyBodyForce(CurrentModel, &
        'Particle Initialization Weight')
    IF( GotWeight ) THEN
      CALL Info('InitializeParticles','Using weight when creating particles',Level=8)
      ALLOCATE( Weight( Mesh % MaxElementNodes ) )
      Weight = 0.0_dp
    END IF


    !------------------------------------------------------------------------
    ! Use a simple bounding box for initializatin
    ! By default a local bounding box is used...
    !-------------------------------------------------------------------------  
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
    
    
    IF( InitMethod == 'box random cubic' .OR. InitMethod == 'box uniform cubic') THEN
      Diam = 2 * GetCReal( Params,'Particle Cell Radius',GotIt)
      IF(.NOT. GotIt ) THEN
        Diam = GetCReal( Params,'Particle Cell Size',GotIt)
      END IF
      IF(.NOT. GotIt ) THEN
        Diam = 2 * GetCReal( Params,'Particle Radius',GotIt)
      END IF
      IF(.NOT. GotIt ) THEN
        CALL Fatal('InitializeParticles','Size of unit cell not given')
      END IF
      
      nx = NINT ( ( MaxCoord(1) - MinCoord(1) ) / Diam )
      ny = NINT( ( MaxCoord(2) - MinCoord(2) ) / Diam )
      IF( dim == 3 ) THEN
        nz = NINT( ( MaxCoord(3) - MinCoord(3) ) / Diam )
      ELSE
        nz = 1
      END IF
    END IF

    AdvVar => VariableGet( Mesh % Variables,'AdvectorData')
    
    !------------------------------------------------------------------------
    ! Now decide on the number of particles.
    !-------------------------------------------------------------------------  
    IF( PRESENT( AppendParticles ) ) THEN
      Offset = Particles % NumberOfParticles
    ELSE
      Offset = 0
    END IF
    
    
    IF( PRESENT( InitParticles ) ) THEN
      NewParticles = InitParticles
    ELSE IF( ASSOCIATED( AdvVar ) ) THEN
      NewParticles = SIZE( AdvVar % Values )
      CALL Info('InitializeParticles','Using pre-existing ParticleData variable to define particles!')
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
              CALL Fatal('InitializeParticles','Could not determine the number of new particles!')
            END IF
          END IF
        END IF
      END IF
    END IF
    
    IF( ParEnv% PEs == 1 ) THEN
      TotParticles = NewParticles
    ELSE
      TotParticles = NINT( ParallelReduction( 1.0_dp * NewParticles ) )
    END IF
    
    IF( TotParticles == 0 ) THEN
      CALL Fatal('InitializeParticles','No Particles to Initialize')
    ELSE
      WRITE( Message,'(A,I0)') 'Number of Particles: ',TotParticles
      CALL Info('InitializeParticles',Message,Level=6)
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
        CALL Fatal('InitializeParticles','Group used inconsistently!')
      END IF
      Particles % Group(Offset+1:LastParticle) = Group
    END IF
      
    Particles % NumberOfParticles = LastParticle
    
    Velocity => Particles % Velocity
    Coordinate => Particles % Coordinate
    
    
    SELECT CASE ( InitMethod ) 
      
    CASE ('nodal ordered')
      CALL Info('InitializeParticles',&
          'Initializing particles evenly among nodes',Level=10)

      Particles % NumberOfParticles = NewParticles
      DO i=1,NewParticles
        k = Offset + i
        j = (nonodes-1)*(i-1)/(NewParticles-1)+1
        IF( GotMask ) j = InvPerm(j)
        Coordinate(k,1) = Mesh % Nodes % x(j)
        Coordinate(k,2) = Mesh % Nodes % y(j)
        IF( dim == 3 ) Coordinate(k,3) = Mesh % Nodes % z(j)
      END DO


      IF( SaveParticleOrigin ) THEN
        DO i=1,Mesh % NumberOfBulkElements
          CurrentElement => Mesh % Elements(i)
          NodeIndexes =>  CurrentElement % NodeIndexes
          n = CurrentElement % TYPE % NumberOfNodes
          DO j=1,n
            k = NodeIndexes(j)
            IF( GotMask ) THEN
              k = MaskPerm(k)
              IF( k == 0 ) CYCLE
            END IF
            Particles % ElementIndex(k) = i
            Particles % NodeIndex(k) = NodeIndexes(j)
          END DO
        END DO
      END IF


    CASE ('elemental random')
      CALL Info('InitializeParticles',&
          'Initializing particles randomly within elements',Level=10)

      n = Mesh % MaxElementNodes 
      ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )

      Particles % NumberOfParticles = NewParticles

      MaxDetJ = 0.0_dp
      MinDetJ = HUGE( MinDetJ )
      MaxWeight = -HUGE( MaxWeight ) 
      MinWeight = HUGE( MinWeight )

      DO i = 1, NoElements

        j = i
        IF( GotMask ) j = InvPerm(j)

        CurrentElement => Mesh % Elements(j)
        NodeIndexes =>  CurrentElement % NodeIndexes
        n = CurrentElement % TYPE % NumberOfNodes

        ! If weight is used see that we have a weight, and that it is positive
        IF( GotWeight ) THEN
          IF( j > Mesh % NumberOfBulkElements ) CYCLE

          body_id = CurrentElement % BodyId
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

        DetJ = ElementSize( CurrentElement, Nodes ) 
        MaxDetJ = MAX( MaxDetJ, DetJ ) 
        MinDetJ = MIN( MinDetJ, DetJ )
      END DO

      WRITE( Message,'(A,ES12.3)') 'Maximum size of elements:',MaxDetJ
      CALL Info('InitializeParticle',Message,Level=8)
      WRITE( Message,'(A,ES12.3)') 'Minimum size of elements:',MinDetJ
      CALL Info('InitializeParticle',Message,Level=8)
      IF( GotWeight ) THEN
        WRITE( Message,'(A,ES12.3)') 'Maximum weight in elements:',MaxWeight
        CALL Info('InitializeParticle',Message,Level=8)
        WRITE( Message,'(A,ES12.3)') 'Minimum weight in elements:',MinWeight
        CALL Info('InitializeParticle',Message,Level=8)
      END IF

      IF( MaxWeight < 0.0 ) THEN
        CALL Info('InitializeParticle','No positive weight!')
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
                
        CurrentElement => Mesh % Elements(j)
        NodeIndexes =>  CurrentElement % NodeIndexes
        n = CurrentElement % TYPE % NumberOfNodes

        Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        IF( CheckForSize ) THEN
          DetJ = ElementSize( CurrentElement, Nodes ) 

          ! The weight could be computed really using the integration point
          ! Here we assumes constant weight within the whole element. 
          IF( GotWeight ) THEN
            IF( j > Mesh % NumberOfBulkElements ) CYCLE

            body_id = CurrentElement % BodyId
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
        Coord = RandomPointInElement( CurrentElement, Nodes )

        i = i + 1
        k = Offset + i
        Coordinate(k,1:dim) = Coord(1:dim)

        ! Only a bulk element may own a particle
        IF( j <= Mesh % NumberOfBulkElements ) THEN
          Particles % ElementIndex(k) = j
        END IF

        IF( i == NewParticles ) EXIT
      END DO
      DEALLOCATE(Nodes % x, Nodes % y, Nodes % z)

    CASE ('elemental ordered')
      CALL Info('InitializeParticles',&
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
        
        CurrentElement => Mesh % Elements(j)
        NodeIndexes =>  CurrentElement % NodeIndexes
        n = CurrentElement % TYPE % NumberOfNodes
        Coord(1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n
        Coord(2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n
        IF( dim == 3 ) Coord(3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n

        Coordinate(k,1:dim) = Coord(1:dim)
        
        ! Only a bulk element may own a particle
        IF( j <= Mesh % NumberOfBulkElements ) THEN
          Particles % ElementIndex(i) = j
        END IF
      END DO

      IF( SaveParticleOrigin ) THEN
        ! For now the initial index is confusingly named always "NodeIndex"
        Particles % NodeIndex = Particles % ElementIndex
      END IF

    CASE ('advector')
      CALL Info('InitializeParticles',&
          'Initializing particles evenly on scaled dg points',Level=10)

      AdvVar => VariableGet( Mesh % Variables,'AdvectorData' )
      IF( .NOT. ASSOCIATED( AdvVar ) ) THEN
        CALL Fatal('InitializeParticles','Variable >AdvectorData< should exist!')
      END IF

      VariableName = ListGetString( Params,'Velocity Variable Name',GotIt )
      IF(.NOT. GotIt) VariableName = 'Flow Solution'
      Var => VariableGet( Mesh % Variables, TRIM(VariableName) )
      IF( .NOT. ASSOCIATED( Var ) ) THEN
        CALL Fatal('InitializeParticles','Velocity variable needed to initialize advector')
      END IF
      vdofs = Var % Dofs
      
      NewParticles = SIZE( AdvVar % Values )
      
      DO i=1,NoElements                       
        CurrentElement => Mesh % Elements(i)
        NodeIndexes =>  CurrentElement % NodeIndexes
        n = CurrentElement % TYPE % NumberOfNodes

        IF( AdvVar % TYPE /= Variable_on_gauss_points ) THEN        
          Center(1) = SUM( Mesh % Nodes % x(NodeIndexes ) ) / n
          Center(2) = SUM( Mesh % Nodes % y(NodeIndexes ) ) / n
          IF( dim == 3 ) Center(3) = SUM( Mesh % Nodes % z(NodeIndexes ) ) / n
          DO j=1,dim
            CenterVelo(j) = SUM( Var % Values(vdofs*(Var % Perm(NodeIndexes)-1)+j) ) / n
          END DO
        END IF

        IF( AdvVar % TYPE == Variable_on_elements ) THEN
          No = AdvVar % Perm( i )
          IF( No == 0 ) CYCLE
          Coordinate(No,1:dim) = Center(1:dim)

          Velocity(No,1:dim) = CenterVelo(1:dim) 

          Particles % ElementIndex(No) = i
          IF( SaveParticleOrigin ) THEN
            Particles % NodeIndex(No) = No
          END IF


        ELSE IF( AdvVar % Type == Variable_on_nodes_on_elements ) THEN

          BLOCK
            REAL(KIND=dp) :: DgScale
            LOGICAL :: GotScale
          
            DGScale = ListGetCReal( Params,'DG Nodes Scale',GotScale )
            IF(.NOT. GotScale ) DgScale = 1.0 / SQRT( 3.0_dp ) 
            GotScale = ( ABS( DGScale - 1.0_dp ) > TINY( DgScale ) )
            
            DO j = 1, n
              No = AdvVar % Perm( CurrentElement % DgIndexes(j) )
              IF( No == 0 ) CYCLE
              k = NodeIndexes(j)
              Coord(1) = Mesh % Nodes % x(k)
              Coord(2) = Mesh % Nodes % y(k)
              IF( dim == 3 ) Coord(3) = Mesh % Nodes % z(k)            

              DO l=1,dim
                Velo(l) = Var % Values(vdofs*(Var % Perm(k)-1)+l) 
              END DO

              IF( GotScale ) THEN
                Coord(1:dim) = Center(1:dim) + ( Coord(1:dim) - Center(1:dim) ) * DgScale 
                Velo(1:dim) = CenterVelo(1:dim) + (Velo(1:dim) - CenterVelo(1:dim)) * DgScale
              END IF

              Coordinate(No,1:dim) = Coord(1:dim) 
              Velocity(No,1:dim) = Velo(1:dim) 

              Particles % ElementIndex(No) = i

              IF( SaveParticleOrigin ) THEN
                Particles % NodeIndex(No) = No                       
              END IF
            END DO
          END BLOCK

            
        ELSE IF( AdvVar % TYPE == Variable_on_gauss_points ) THEN

          BLOCK
            TYPE(GaussIntegrationPoints_t) :: IP
            REAL(KIND=dp) :: detJ, Basis(27)
            LOGICAL :: stat
            TYPE(Nodes_t), SAVE :: Nodes
            INTEGER :: m
            LOGICAL :: Debug

            Debug = ( i == 0 )
            
            IF( i == 1 ) THEN
              m = 27
              ALLOCATE( Nodes % x(m), Nodes % y(m), Nodes % z(m))
            END IF

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

            IP = GaussPoints(CurrentElement, m )
            
            DO j = 1, IP % n 
              stat = ElementInfo( CurrentElement, Nodes, IP % v(j), IP % u(j), IP % w(j), detJ, Basis ) 
              No = AdvVar % Perm(i) + j
              
              Coord(1) = SUM( Basis(1:n) * Nodes % x(1:n) )
              Coord(2) = SUM( Basis(1:n) * Nodes % y(1:n) )
              IF( dim == 3 ) Coord(3) = SUM( Basis(1:n) * Nodes % z(1:n) )

              DO l=1,dim
                Velo(l) = SUM( Basis(1:n) * Var % Values(vdofs*(Var % Perm(NodeIndexes)-1)+l ) )
              END DO

              IF( Debug ) THEN
                PRINT *,'j:',j,m,No,Coord(1:dim),Velo(1:dim)
              END IF
                
              Coordinate(No,1:dim) = Coord(1:dim)
              Velocity(No,1:dim) = Velo(1:dim)
              
              Particles % ElementIndex(No) = i
              IF( SaveParticleOrigin ) THEN
                Particles % NodeIndex(No) = No                       
              END IF
            END DO            
          END BLOCK

        END IF
      END DO
      
      
    CASE ('sphere random')
      CALL Info('InitializeParticles',&
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
      CALL Info('InitializeParticles',&
          'Initializing particles randomly in a grid',Level=10)

      nmax = nx * ny * nz
      IF( nmax < NewParticles ) THEN
        CALL Fatal('InitializeParticles','More particles than places in unit cell')
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
      CALL Info('InitializeParticles',&
          'Initializing particles in a grid',Level=10)

      nmax = nx * ny * nz
      IF( nmax /= NewParticles ) THEN
        CALL Fatal('InitializeParticles','Wrong number of particles')
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
      CALL Info('InitializeParticles',&
          'Initializing particles using given coordinates',Level=10)

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


    IF( GotMask .AND. ASSOCIATED(InvPerm) ) DEALLOCATE( InvPerm ) 

    !------------------------------------------------------------------------
    ! Velocities may be initialized using a given list, or obtaining them
    ! from random even or maxwell boltzmann distributions. These are additive to 
    ! allow bulk velocities with the random one.
    !-------------------------------------------------------------------------

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
    coeff = ListGetCReal( Params,'Initial Velocity Amplitude',GotIt)


    SELECT CASE ( InitMethod ) 

    CASE ('nodal velocity')
      CALL Info('InitializeParticles',&
          'Initializing velocities from the corresponding nodal velocity',Level=10)
      
      VariableName = ListGetString( Params,'Velocity Variable Name',GotIt )
      IF(.NOT. GotIt) VariableName = 'Flow Solution'
      Var => VariableGet( Mesh % Variables, TRIM(VariableName) )
      IF( .NOT. ASSOCIATED( Var ) ) THEN
        CALL Fatal('InitializeParticles','Velocity variable needed for method >nodal velocity<')
      END IF

      vdofs = Var % Dofs
      DO i=1,NewParticles
        k = Offset + i
        l = Particles % NodeIndex(i)
        l = Var % Perm(l)
        DO j=1,dim
          Velocity(k,j) = Var % Values(vdofs*(l-1)+j)
        END DO
      END DO

    CASE ('elemental velocity') 
      CALL Info('InitializeParticles',&
          'Initializing velocities from the corresponding elemental velocity',Level=10)
      
      VariableName = ListGetString( Params,'Velocity Variable Name',GotIt )
      IF(.NOT. GotIt) VariableName = 'Flow Solution'
      Var => VariableGet( Mesh % Variables, TRIM(VariableName) )
      IF( .NOT. ASSOCIATED( Var ) ) THEN
        CALL Fatal('InitializeParticles','Velocity variable needed for method >elemental velocity<')
      END IF
      vdofs = Var % Dofs

      DO i=1,NewParticles
        k = Offset + i
        l = Particles % NodeIndex(i) ! now an elemental index
        NodeIndexes => Mesh % Elements(l) % NodeIndexes
        DO j=1,dim
          Velocity(k,j) = SUM( Var % Values(vdofs*(Var % Perm(NodeIndexes)-1)+j) ) / SIZE( NodeIndexes )
        END DO
      END DO

    CASE ('advector') 
      CALL Info('InitializeParticles',&
          'Velocities have been initialized together with the position',Level=10)
      
    CASE ('thermal random')  
       CALL Info('InitializeParticles',&
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
      CALL Info('InitializeParticles',&
          'Initializing velocities from a even distribution',Level=10)
      
      DO i=1,NewParticles
        k = Offset + i
        DO j=1,dim
          Velo(j) = (2*EvenRandom()-1)
        END DO
        Velocity(k,:) = Velocity(k,:) + coeff * Velo(1:dim)
      END DO

    CASE ('constant random')
      CALL Info('InitializeParticles',&
          'Initializing constant velocities with random direction',Level=10)
      
      DO i=1,NewParticles
        k = Offset + i
        DO j=1,dim
          Velo(j) =  (2*EvenRandom()-1)
        END DO
        Velo(1:dim) = Velo(1:dim) / SQRT(SUM(Velo(1:dim)**2))
        Velocity(k,:) = Velocity(k,:) + coeff * Velo(1:dim)
      END DO

    CASE ('constant 2d')
      CALL Info('InitializeParticles',&
          'Initializing constant velocities evenly to space',Level=10)
      
      DO i=1,NewParticles
        k = Offset + i
        Phi = 2.0_dp * PI * i / NewParticles
        Velo(1) = coeff * COS( Phi )
        Velo(2) = coeff * SIN( Phi )
        Velocity(k,:) = Velocity(k,:) + Velo(1:dim)
      END DO

    CASE DEFAULT

    END SELECT


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
    END IF
    IF( ListCheckPresentAnyBodyForce( CurrentModel,&
        'Particle Time Integral Source') ) THEN
      CALL ParticleVariableCreate( Particles,'Particle Time Integral' )
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
  !> consistant strategy with the above algorithm.
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
    INTEGER :: i,j,k,n,FaceIndex,MaxTrials,bc_id,cons_id,ElementIndex0,ParticleStatus0, &
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
        DebugNo, DebugPart
    
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

      Eps = ListGetConstReal( Params,'Particle Hit Tolerance',Stat)
      IF(.NOT. Stat) Eps = 1.0e-10
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
        WRITE(Message,'(A,3ES10.3)') 'Losing particle '//TRIM(I2S(No))//' in: ',Rfin(1:3)
        CALL Info('LocateParticlesInMesh',Message,Level=15)
        
        ParticleStatus = PARTICLE_LOST
        EXIT
      END IF

      PrevElement => Element
      Element => NextElement
    END DO

    IF( i >= MaxTrials ) THEN
      PRINT *,'Used maximum number of trials',MaxTrials,No
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
    REAL(KIND=dp) :: Rinit(3), Rfin(3),Rfin0(3),Velo(3), Velo0(3), dtime
    LOGICAL :: Stat, InitLocation, AccurateAtFace, AccurateAlways, AccurateNow, debug
    INTEGER :: FaceIndex, FaceIndex0, Status0, InitStatus
    REAL(KIND=dp) :: Lambda
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: DtVar
    
    CALL Info('LocateParticles','Locating particles in mesh',Level=10)

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
        CALL Fatal('ParticleAdvanceTimesteo','Variable timestep, > particle dt < should exist!')
      END IF
    END IF


100 NoParticles = Particles % NumberOfParticles

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

200   ElementIndex = GetParticleElement( Particles, No )
      Rfin = GetParticleCoord( Particles, No )
      Velo = GetParticleVelo( Particles, No )      
      IF( AccurateNow ) Rinit = GetParticlePrevCoord( Particles, No )        
      Rinit = GetParticlePrevCoord( Particles, No )        

      IF( debug ) THEN
        PRINT *,parenv % mype, 'going No',No,'Element',ElementIndex,'Face',FaceIndex,'Status',Status
        PRINT *,parenv % mype, 'going Init:    ',Rinit(1:dim),Rfin(1:dim)
        PRINT *,parenv % mype, 'going Velo:',GetParticleVelo(Particles,No), Velo(1:dim)
      END IF
      
      CALL LocateParticleInMeshMarch(ElementIndex, Rinit, Rfin, InitLocation, &
          Status,AccurateNow, FaceIndex, Lambda, Velo, No, ParticleWallKernel, Particles )

      IF( debug ) THEN
        PRINT *,parenv % mype, 'leving No',No,'Element',ElementIndex,'Face',FaceIndex,'Status',Status
        PRINT *,parenv % mype, 'leving Init:    ',Rinit(1:dim),Rfin(1:dim)
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
          CALL Warn('LocateParticles','Difference between robust and accurate?')
        END IF   

        Status = Status0
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

    ! Change the partion in where the particles are located
    ! Only applies to parallel cases.
    !------------------------------------------------------------------------
    PartitionChanges = ChangeParticlePartition( Particles )
    IF( PartitionChanges > 0 ) THEN
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
  FUNCTION ParticleElementInfo( CurrentElement, GlobalCoord, &
      SqrtElementMetric, Basis, dBasisdx ) RESULT ( stat )
    
    TYPE(Element_t), POINTER :: CurrentElement
    REAL(KIND=dp) :: GlobalCoord(:), SqrtElementMetric, LocalDistance
    REAL(KIND=dp) :: Basis(:)
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:)
    LOGICAL :: Stat, Debug
    INTEGER :: Misses(2) = 0
  
    SAVE Misses    


    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: LocalCoord(3),u,v,w
    INTEGER :: n
    
    SAVE ElementNodes
    
    n = CurrentElement % TYPE % NumberOfNodes
    CALL GetElementNodes(ElementNodes,CurrentElement)
    
    Stat = PointInElement( CurrentElement, ElementNodes, &
        GlobalCoord, LocalCoord, GlobalEps = -1.0_dp, LocalEps = 1.0e3_dp, &
	LocalDistance = LocalDistance ) 

    IF( .NOT. Stat ) THEN
      Misses(1) = Misses(1) + 1

      IF( MODULO( SUM( Misses ), 101 ) == 100 ) PRINT *,'Misses:',Misses

      IF( .FALSE.) THEN
        IF( .NOT. Stat ) THEN
          CALL Warn('ParticleElementInfo','Should have found the node!')
        ELSE
          CALL Warn('ParticleElementInfo','Distance from element higher than expected!')
        END IF
        PRINT *,'LocalDistance:',LocalDistance,'Element:',CurrentElement % ElementIndex
        PRINT *,'Nodes X:',ElementNodes % x(1:n) - GlobalCoord(1)
        PRINT *,'Nodes Y:',ElementNodes % y(1:n) - GlobalCoord(2)
        PRINT *,'Nodes Z:',ElementNodes % z(1:n) - GlobalCoord(3)
      END IF
      RETURN
    END IF

    u = LocalCoord(1)
    v = LocalCoord(2)
    w = LocalCoord(3)
    
    stat = ElementInfo( CurrentElement, ElementNodes, U, V, W, SqrtElementMetric, &
        Basis, dBasisdx )
    IF(.NOT. Stat) Misses(2) = Misses(2) + 1
    
  END FUNCTION ParticleElementInfo
  


  !-------------------------------------------------------------------------
  !> The routine returns velocity and optionally a gradient of velocity.
  !> These kind of functions are needed repeated and therefore to reduced the 
  !> size of individual solvers it has been hard coded here. 
  !--------------------------------------------------------------------------
  
  SUBROUTINE GetVectorFieldInMesh(Var, CurrentElement, Basis, Velo, dBasisdx, GradVelo )
    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t) :: CurrentElement
    REAL(KIND=dp) :: Basis(:), Velo(:) 
    REAL(KIND=dp), OPTIONAL :: dBasisdx(:,:), GradVelo(:,:)
    
    TYPE(Valuelist_t), POINTER :: Params
    INTEGER, POINTER :: LocalPerm(:)
    REAL(KIND=dp), POINTER :: LocalVelo(:,:)
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: VeloFieldDofs
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
    
    n = CurrentElement % TYPE % NumberOfNodes
    LocalPerm(1:n) = Var % Perm( CurrentElement % NodeIndexes )
    npos = COUNT ( LocalPerm(1:n) > 0 )
    
    
    IF( npos == 0 ) RETURN
    
    !-----------------------------------------------------------------
    ! compute the velocity also for case when the particle
    ! has just crossed the boundary. For example, its floating on the 
    ! fluid boundary. This is a little bit fishy and could perhaps 
    ! only be done conditionally....
    ! Can't really determine the gradient here
    !-----------------------------------------------------------------
    VeloFieldDofs = Var % Dofs
    IF( npos == n ) THEN
      DO i=1,n
        j = LocalPerm(i)
	DO k=1,dim
          LocalVelo(i,k) = Var % Values( VeloFieldDofs*(j-1)+k)
        END DO
      END DO
    ELSE    
      IF(.NOT. InterfaceNodes ) RETURN

      SumBasis = 0.0_dp
      DO i=1,n
        j = LocalPerm(i)
        IF( j > 0 ) THEN
          SumBasis = SumBasis + Basis(i)
          DO k=1,dim
            LocalVelo(i,k) = Var % Values( VeloFieldDofs*(j-1)+k)
          END DO
        ELSE
          Basis(i) = 0.0_dp
          LocalVelo(i,1:dim) = 0.0_dp
        END IF
      END DO
    END IF
    

    DO i=1,dim
      Velo(i) = SUM( Basis(1:n) * LocalVelo(1:n,i) )
      IF( PRESENT( GradVelo ) ) THEN
        DO j=1,dim
          GradVelo(i,j) = SUM( dBasisdx(1:n,j) * LocalVelo(1:n,i) )
        END DO
      END IF
    END DO
    
    IF( npos < n ) THEN
      Velo(1:dim) = Velo(1:dim) / SumBasis
      IF( PRESENT( GradVelo ) ) THEN
        GradVelo(:,1:dim) = GradVelo(:,1:dim) / SumBasis
      END IF
    END IF
   
  END SUBROUTINE GetVectorFieldInMesh


  !-------------------------------------------------------------------------
  !> The routine returns a potential and its gradient.
  !--------------------------------------------------------------------------
  
  SUBROUTINE GetScalarFieldInMesh(Var, CurrentElement, Basis, Pot, dBasisdx, GradPot )
    
    TYPE(Variable_t), POINTER :: Var
    TYPE(Element_t) :: CurrentElement
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
    
    n = CurrentElement % TYPE % NumberOfNodes
    IF( ASSOCIATED( Var % Perm ) ) THEN
      LocalPerm(1:n) = Var % Perm( CurrentElement % NodeIndexes )
      IF( .NOT. ALL ( LocalPerm(1:n) > 0 )) RETURN
      LocalField(1:n) = Var % Values( LocalPerm(1:n) )
    ELSE
      ! Some variables do not have permutation, most importantly the node coordinates
      LocalField(1:n) = Var % Values( CurrentElement % NodeIndexes )
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
  !> Also if the particle is splitted between two elements then this 
  !> routine can assess the data on the secondary mesh.
  !-------------------------------------------------------------
  FUNCTION GetMaterialPropertyInMesh(PropertyName, BulkElement, Basis, &
      BulkElement2, VolumeFraction ) RESULT ( Property )
    
    CHARACTER(LEN=MAX_NAME_LEN) :: PropertyName
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
    
    ! And finally, add the closest neigbours to the table 
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
      IF ( .NOT. PI % INTERFACE(node) ) CYCLE
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
      IF ( .NOT. PI % INTERFACE(node) ) CYCLE
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
            CALL Info('GetNextNeighbour','Allocating more space: '//TRIM(I2S(ListSize)))
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
  SUBROUTINE ParticleAdvanceTimestep( Particles, RKstep )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: RKstep

    REAL(KIND=dp) :: dtime
    TYPE(Variable_t), POINTER :: Var, TimeVar, DistVar, DtVar
    LOGICAL :: GotVar, GotTimeVar, GotDistVar, MovingMesh
    REAL(KIND=dp) :: ds, dCoord(3),Coord(3),Velo(3),Speed0,Speed
    INTEGER :: dim, Status, TimeOrder, No, NoMoving
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoParticles
    LOGICAL :: Found, Visited = .FALSE.,RK2,HaveSpeed0

    REAL(KIND=dp) :: mass, drag
    REAL(KIND=dp), POINTER :: massv(:), dragv(:)
    LOGICAL :: GotMass, GotDrag
    INTEGER :: CurrGroup, PrevGroup, NoGroups
    
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
            CALL Fatal('ParticleAdvanceTimesteo','Variable timestep, > particle dt < should exist!')
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
      IF(.NOT. GotMass) CALL Fatal('ParticleAdvanceTime',&
          '> Particle Mass < should be given!')
    ELSE IF( TimeOrder == 1 ) THEN
      IF( NoGroups > 1 ) THEN
        dragv => ListGetConstRealArray1( Params,'Particle Drag Coefficient',GotDrag)
        Drag = 0.0_dp
      ELSE
        Drag = ListGetConstReal( Params,'Particle Drag Coefficient',GotDrag)
      END IF
      IF(.NOT. GotDrag) CALL Fatal('ParticleAdvanceTime',&
          '> Particle Drag Coefficient < should be given!')
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
        CALL Fatal('ParticleAdvanceTimestep','Unknown time order')
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
  SUBROUTINE ParticlePathIntegral( Particles, RKstep )
    TYPE(Particle_t), POINTER :: Particles
    INTEGER, OPTIONAL :: RKstep

    TYPE(Variable_t), POINTER :: TimeIntegVar, DistIntegVar, DtVar
    LOGICAL :: GotVar, RK2
    REAL(KIND=dp) :: ds,dtime,Coord(3),PrevCoord(3),LocalCoord(3),Velo(3),u,v,w,&
        SourceAtPath,detJ,RKCoeff
    INTEGER :: dim, Status
    TYPE(ValueList_t), POINTER :: Params
    INTEGER :: NoParticles, No, n, NoVar, i, j, bf_id
    LOGICAL :: Found, Stat, Visited = .FALSE.
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), POINTER :: Basis(:), Source(:), dBasisdx(:,:)
    INTEGER, POINTER :: Indexes(:)
    TYPE(ValueList_t), POINTER :: BodyForce
    TYPE(Element_t), POINTER :: Element
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: str, VariableName
    LOGICAL :: TimeInteg, DistInteg, UseGradSource


    SAVE TimeInteg, DistInteg, dim, Visited, Mesh, DtVar, Basis, Source, Nodes, Params, &
        TimeIntegVar, DistIntegVar, UseGradSource, dBasisdx

    CALL Info('ParticlePathIntegral','Integrating variables over the path',Level=12)


    ! If Runge-Kutta is used take the mid-point rule.
    IF( RKStep > 1 ) RETURN
    
    IF(.NOT. Visited ) THEN
      Visited = .TRUE.

      Params => ListGetSolverParams()
      Mesh => CurrentModel % Solver % Mesh
      dim = Particles % dim

      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), Source(n), Nodes % x(n), Nodes % y(n), Nodes % z(n), &
          dBasisdx(n,3) )
      Basis = 0.0_dp
      Source = 0.0_dp
      Nodes % x = 0.0_dp
      Nodes % y = 0.0_dp
      Nodes % z = 0.0_dp
      dBasisdx = 0.0_dp

      UseGradSource = GetLogical( Params,'Source Gradient Correction',Found)
      ! If the correction is not given follow the logic of velocity estimation
      IF(UseGradSource .AND. Particles % RK2 ) THEN
        CALL Warn('ParticlePathIntegral','Quadratic source correction incompatibe with Runge-Kutta')
        UseGradSource = .FALSE.
      END IF

      IF( .NOT. Particles % DtConstant ) THEN
        DtVar => ParticleVariableGet( Particles,'particle dt')
        IF(.NOT. ASSOCIATED( DtVar ) ) THEN
          CALL Fatal('ParticleAdvanceTimesteo','Variable timestep, > particle dt < should exist!')
        END IF
      END IF
      
      TimeIntegVar => ParticleVariableGet( Particles,'particle time integral')
      TimeInteg = ASSOCIATED( TimeIntegVar )

      DistIntegVar => ParticleVariableGet( Particles,'particle distance integral')
      DistInteg = ASSOCIATED( DistIntegVar )
    END IF

    ! Nothing to integrate over
    IF( .NOT. (TimeInteg .OR. DistInteg ) ) RETURN

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

      Coord = 0._dp
      Coord(1:dim) = Particles % Coordinate(No,:) 
      Velo  = 0._dp
      Velo(1:dim) = Particles % Velocity(No,:) 	
      
      Element => Mesh % Elements( Particles % ElementIndex(No) )            
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
        Source(1:n) = ListGetReal( BodyForce,'Particle Time Integral Source', &
            n, Indexes, Found )
        IF( Found ) THEN
          SourceAtPath = SUM( Basis(1:n) * Source(1:n) )
          IF( UseGradSource ) THEN
            DO i=1,dim
              SourceAtPath = SourceAtPath + 0.5*SUM( dBasisdx(1:n,i) * Source(1:n) ) * &
                  ( PrevCoord(i) - Coord(i) )
            END DO
          END IF
          TimeIntegVar % Values(No) = TimeIntegVar % Values(No) + dtime * SourceAtPath
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
      
        Source(1:n) = ListGetReal( BodyForce,'Particle Distance Integral Source', &
            n, Indexes, Found ) 
        IF( Found ) THEN
          SourceAtPath = SUM( Basis(1:n) * Source(1:n) )
          IF( UseGradSource ) THEN
            DO i=1,dim
              SourceAtPath = SourceAtPath + 0.5*SUM( dBasisdx(1:n,i) * Source(1:n) ) * &
                  ( PrevCoord(i) - Coord(i) )
            END DO
          END IF
          DistIntegVar % Values(No) = DistIntegVar % Values(No) + ds * SourceAtPath
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
    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    
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
!> Set a the timestep for the particles.
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
        CharSpeed, CharTime, dtave
    LOGICAL :: GotIt,TfinIs,NStepIs,DsGoalIs,HgoalIs,DtIs
    INTEGER :: nstep, TimeStep, PrevTimeStep = -1, flag
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: TimeVar, DtVar
    
    SAVE dt0,dsgoal,hgoal,dtmax,dtmin,DtIs,Nstep,&
        tprev,Tfin,TfinIs,DsGoalIs,HgoalIs,PrevTimeStep, &
	DtVar,TimeVar
    
    dtout = 0.0_dp

    IF( InitInterval ) THEN
      Params => ListGetSolverParams()
      
      ! directly defined timestep
      dt0 = GetCReal(Params,'Timestep Size',DtIs)
      
      ! Constraint by absolute step size taken (in length units)
      dsgoal = GetCReal( Params,'Timestep Distance',DsGoalIs)
      
      ! Constraint by relative step size taken (1 means size of the element)
      hgoal = GetCReal( Params,'Timestep Courant Number',HGoalIs)
      
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

        !PRINT *,'ratio of timesteps:',tfin/dt
      ELSE IF( tfinIs ) THEN
        dt = tfin / Nstep
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
      tprev = tprev + dt
      Particles % dtime = dt
      dtout = dt
    ELSE 
      DtVar % Values = 0.0_dp
      dtave = 0.0_dp
      
      DO No = 1, Particles % NumberOfParticles

        Status = Particles % Status( No )
        IF ( Status >= PARTICLE_LOST ) CYCLE
        IF ( Status <= PARTICLE_INITIATED ) CYCLE
        IF ( Status == PARTICLE_WALLBOUNDARY ) CYCLE
        IF ( Status == PARTICLE_FIXEDCOORD ) CYCLE

	tprev = TimeVar % Values(No)        

	flag = 1

        IF( DtIs ) THEN
          dt = dt0 
	  flag = 2
        ELSE IF( DsGoalIs ) THEN
          CharSpeed = CharacteristicSpeed( Particles, No )     
          dt = dsgoal / CharSpeed
	  flag = 3
        ELSE IF( HgoalIs ) THEN
          CharTime = CharacteristicElementTime( Particles, No )     
          dt = Hgoal * CharTime ! ElementH / Speed
	  flag = 4
        ELSE IF( tfinIs ) THEN
          dt = tfin / Nstep
	  flag = 5
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

        dtout = MAX( dtout, dt )
        dtave = dtave + dt
      END DO

      dtave = dtave / Particles % NumberOfParticles

      WRITE(Message,'(A,ES12.3)') 'Average particle timestep:',dtave
      CALL Info('GetParticleTimestep', Message,Level=12)           
    END IF
    
    dtout = ParallelReduction( dtout, 2 )       
      
    WRITE(Message,'(A,ES12.3)') 'Maximum particle timestep:',dtout
    CALL Info('GetParticleTimestep', Message,Level=12)           
  
    
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
!>  Given a variable name, get a handle to it.
!------------------------------------------------------------------------------
  FUNCTION ParticleVariableGet( Particles, Name ) RESULT ( Var )
    
    TYPE(Particle_t), POINTER :: Particles
    CHARACTER(LEN=*) :: Name
    TYPE(Variable_t), POINTER :: Var
!------------------------------------------------------------------------------
    Var => VariableGet( Particles % Variables, Name )

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
        CALL Info('ParticleVariablesResize','Using compact size of: '//TRIM(I2S(oldsize)),Level=12)
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
    INTEGER :: i,j,n,dofs,Vari,Rank,dim, NoParticles, MinSaveStatus, MaxSaveStatus
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


    IF( VisitedTimes == 1 ) THEN
      Params => ListGetSolverParams()
      FilePrefix = ListGetString(Params,'Filename Prefix')
      CALL WriteParticleFileNames(FilePrefix, dim)
      
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

    CALL Info('ParticleOutputTable','Saving at maximum '//TRIM(I2S(NoParticles))//' particles',Level=6)
    
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
      
      CHARACTER(LEN=MAX_NAME_LEN) :: Prefix
      INTEGER :: dim

      CHARACTER(LEN=MAX_NAME_LEN) :: FileName
      INTEGER :: i,j,dofs
      TYPE(Variable_t), POINTER :: Solution
      LOGICAL :: ComponentVector, ThisOnly = .TRUE.
      CHARACTER(LEN=1024) :: Txt, FieldName


      WRITE( FileName,'(A,A)') TRIM(FilePrefix),'.dat.names'
      
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
              WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
            ELSE
              WRITE(Txt,'(A,I0)') 'Vector Field ',Vari             
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
                  WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                  CALL Warn('WriteParticleLine', Txt)
                  CYCLE
                END IF
              END IF
            ELSE
              IF(.NOT. ASSOCIATED(Solution)) THEN
                WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                CALL Warn('WriteParticleLine', Txt)
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
               WRITE( TableUnit, '(I2,A)' )  i+j,': '//TRIM(FieldName)//'_'//TRIM(I2S(j))  
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
      
      CHARACTER(LEN=MAX_NAME_LEN) :: Prefix
      INTEGER :: FileNo
      LOGICAL, SAVE :: Visited = .FALSE.


      CHARACTER(LEN=MAX_NAME_LEN) :: FileName
      
      IF( FileNo == 0 ) THEN
        WRITE( FileName,'(A,A)') TRIM(FilePrefix),'.dat'
        IF( .NOT. Visited ) THEN
          CALL Info( 'ParticleOutputTable', 'Saving particle data to file: '//TRIM(FileName), Level=4 )
        END IF
      ELSE
        IF ( FileNo==1 .AND.  .NOT. Visited ) THEN
          WRITE( Message, * ) 'Saving particle data to files: ', TRIM(FilePrefix)//'_*.dat'
          CALL Info( 'ParticleOutputTable', Message, Level=4 )
        END IF
        FileName=TRIM(FilePrefix)//'_'//TRIM(i2s(fileno))//'.dat'
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
      CHARACTER(LEN=1024) :: Txt, FieldName
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
              WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
            ELSE
              WRITE(Txt,'(A,I0)') 'Vector Field ',Vari             
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
                  WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                  CALL Warn('WriteParticleLine', Txt)
                  CYCLE
                END IF
              END IF
            ELSE
              IF(.NOT. ASSOCIATED(Solution)) THEN
                WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                CALL Warn('WriteParticleLine', Txt)
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
    CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix, FileNameGmsh, FileNameOut
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

    WRITE( FileNameGmsh,'(A,A)') TRIM(FilePrefix),'.pos'
    
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
    
    CHARACTER(MAX_NAME_LEN), SAVE :: FilePrefix
    CHARACTER(MAX_NAME_LEN) :: VtuFile, PvtuFile 
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i, j, k, Partitions, Part, ExtCount, FileindexOffSet, &
        Status, MinSaveStatus, MaxSaveStatus, PrecBits, PrecSize, IntSize, &
        iTime 
    CHARACTER(MAX_NAME_LEN) :: Dir
    REAL(KIND=dp) :: SaveNodeFraction, LocalVal(3)
    LOGICAL :: BinaryOutput,AsciiOutput,Found,Visited = .FALSE.,SaveFields
    REAL(KIND=dp) :: DoubleWrk
    REAL :: SingleWrk

    CHARACTER(MAX_NAME_LEN) :: Str
    INTEGER :: NumberOfNodes, ParallelNodes, Dim
    
    SAVE :: MinSaveStatus, MaxSaveStatus
    
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
      CALL Info('VtuOutputSolver','Using single precision arithmetics in output!',Level=7)
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

    NumberOfNodes = 0
    DO i=1,Particles % NumberOfParticles
      IF ( Particles % Status(i) > MaxSaveStatus .OR. &
          Particles % Status(i) < MinSaveStatus )  CYCLE
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


    IF (LEN_TRIM(Mesh % Name) > 0 ) THEN
      Dir = TRIM(Mesh % Name) // "/"
    ELSE
      Dir = "./"
    END IF
    
    IF(Parallel .AND. Part == 0) THEN
      IF( iTime < 10000 ) THEN
        WRITE( PvtuFile,'(A,A,I4.4,".pvtu")' ) TRIM(Dir),TRIM(FilePrefix),iTime
      ELSE
        WRITE( PvtuFile,'(A,A,I0,".pvtu")' ) TRIM(Dir),TRIM(FilePrefix),iTime
      END IF
      CALL WritePvtuFile( PvtuFile )
    END IF
    
    IF ( Parallel ) THEN
      IF( iTime < 10000 ) THEN
        WRITE( VtuFile,'(A,A,I4.4,A,I4.4,".vtu")' ) TRIM(Dir),TRIM(FilePrefix),Part+1,"par",&
            iTime
      ELSE
        WRITE( VtuFile,'(A,A,I4.4,A,I0,".vtu")' ) TRIM(Dir),TRIM(FilePrefix),Part+1,"par",&
            iTime
      END IF
    ELSE
      IF( iTime < 10000 ) THEN
        WRITE( VtuFile,'(A,A,I4.4,".vtu")' ) TRIM(Dir),TRIM(FilePrefix),iTime
      ELSE
        WRITE( VtuFile,'(A,A,I0,".vtu")' ) TRIM(Dir),TRIM(FilePrefix),iTime
      END IF
    END IF

    CALL Info('ParticleOutputVtu','Saving particles to file: '//TRIM(VtuFile),Level=8)
    CALL WriteVtuFile( VtuFile )
    

  CONTAINS

    
    SUBROUTINE WriteVtuFile( VtuFile )
      CHARACTER(LEN=*), INTENT(IN) :: VtuFile
      INTEGER, PARAMETER :: VtuUnit = 58
      TYPE(Variable_t), POINTER :: Var, Solution
      CHARACTER(LEN=512) :: str
      INTEGER :: i,j,k,dofs,Rank,cumn,n,vari,sdofs,IsVector,Offset,PartDim
      CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
          FieldName2, BaseString, OutStr
      CHARACTER :: lf
      LOGICAL :: ScalarsExist, VectorsExist, Found, ParticleMode, ComponentVector, &
          ComplementExists, ThisOnly, Stat
      LOGICAL :: WriteData, WriteXML, Buffered, IsDG   
      INTEGER, POINTER :: Perm(:), Perm2(:), Indexes(:)
      INTEGER, ALLOCATABLE :: ElemInd(:),ElemInd2(:)
      REAL(KIND=dp), POINTER :: Values(:),Values2(:),&
          Values3(:),VecValues(:,:),Basis(:)
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

      IF( SaveFields ) THEN
        
        DO IsVector = 0, 1

          DO Vari = 1, 999

            IF( IsVector == 0 ) THEN
              BaseString = 'Scalar Field'
            ELSE
              BaseString = 'Vector Field'
            END IF

            WRITE(Txt,'(A)') TRIM(BaseString)//' '//TRIM(I2S(Vari))
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
                WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                CALL Warn('WriteVtuXMLFile', Txt)
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
              WRITE(Txt,'(A,I0,A)') TRIM(BaseString)//' ',Vari,' Complement'

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
              IF( IsVector == 1) THEN
                dofs = PartDim
                IF( FieldName == 'velocity' ) THEN
                  VecValues => Particles % Velocity
                ELSE IF( FieldName == 'force') THEN
                  VecValues => Particles % Force 
                ELSE
                  WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                  CALL Warn('WriteVtuXMLFile', Txt)
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
                ELSE
                  WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
                  CALL Warn('WriteVtuXMLFile', Txt)
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
                j = j + 1

                LocalVal = 0.0_dp

                IF( ParticleMode ) THEN
                  IF( IsVector == 1) THEN
                    dofs = dim
                    LocalVal(1:dofs) = VecValues(i,1:dim)
                  ELSE
                    dofs = 1
                    LocalVal(1) = Values(i)
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

          !PRINT *,'ParticleMode:',ParticleMode,j,dim,LocalVal
          
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
      CHARACTER(LEN=512) :: str
      INTEGER :: i,j,k,dofs,Rank,cumn,n,vari,sdofs
      CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
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
        WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
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
            WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
            CALL Warn('WriteVtuXMLFile', Txt)
            CYCLE
          END IF
        END IF
        
        IF( AsciiOutput ) THEN
          WRITE( VtuUnit,'(A)') '      <PDataArray type="Float'//TRIM(I2S(PrecBits))//&
              '" Name="'//TRIM(FieldName)//'" NumberOfComponents="1" format="ascii"/>'    
        ELSE
          WRITE( VtuUnit,'(A)') '      <PDataArray type="Float'//TRIM(I2S(PrecBits))//&
              '" Name="'//TRIM(FieldName)//'" NumberOfComponents="1" format="appended"/>'    
        END IF

      END DO


      !-------------------------------------------------------      
      ! Do the vectors
      !-------------------------------------------------------          

      DO Vari = 1, 99
        WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
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
              WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
              CALL Warn('WriteVtuXMLFile', Txt)
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
            WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
            CALL Warn('WriteVtuXMLFile', Txt)
            CYCLE
          END IF
        END IF

        sdofs = dofs
        IF( AsciiOutput ) THEN
          WRITE( VtuUnit,'(A,I1,A)') '      <PDataArray type="Float'//TRIM(I2S(PrecBits))//'" Name="&
              '//TRIM(FieldName)//'" NumberOfComponents="',sdofs,'" format="ascii"/>'    
        ELSE
          WRITE( VtuUnit,'(A,I1,A)') '      <PDataArray type="Float'//TRIM(I2S(PrecBits))//'" Name="&
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
        WRITE( VtuUnit,'(A)') '      <DataArray type="Float'//TRIM(I2S(PrecBits))//&
            '" NumberOfComponents="3" format="ascii"/>'    
      ELSE
        WRITE( VtuUnit,'(A)') '      <DataArray type="Float'//TRIM(I2S(PrecBits))//&
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

!    USE DefUtils 
!    USE MeshUtils
!    USE ElementDescription
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
    CHARACTER(MAX_NAME_LEN) :: Dir
    REAL(KIND=dp) :: SaveNodeFraction
    LOGICAL :: Found,BinaryOutput,AsciiOutput,SinglePrec,NoFileIndex, &
        Visited = .FALSE.
  
    CHARACTER(MAX_NAME_LEN) :: Str
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
      TYPE(Variable_t), POINTER :: Var, Solution
      CHARACTER(LEN=512) :: str
      INTEGER :: i,j,k,l,dofs,Rank,cumn,n,vari,sdofs,ind,IsVector,IsAppend,GridPoints,Offset
      CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
          FieldName2, BaseString, OutStr
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
          
          WRITE(Txt,'(A)') TRIM(BaseString)//' '//TRIM(I2S(Vari))          
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
            WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
            CALL Warn('WriteVtiXMLFile', Txt)
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
          WRITE(Txt,'(A,I0,A)') TRIM(BaseString),Vari,' Complement'

          FieldName2 = ListGetString( Params, TRIM(Txt), Found )
          IF( Found ) THEN
            Solution => VariableGet( Mesh % Variables, &
                TRIM(FieldName2), ThisOnly )
            IF( ASSOCIATED(Solution)) THEN 
              Values2 => Solution % Values
              Perm2 => Solution % Perm 
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
      CHARACTER(LEN=512) :: str
      INTEGER :: i,j,k,dofs,Rank,cumn,n,vari,sdofs
      CHARACTER(LEN=1024) :: Txt, ScalarFieldName, VectorFieldName, FieldName, &
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
        WRITE(Txt,'(A,I0)') 'Scalar Field ',Vari
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
        WRITE(Txt,'(A,I0)') 'Vector Field ',Vari
        FieldName = ListGetString( Params, TRIM(Txt), Found )
        IF(.NOT. Found) EXIT
        
        Solution => VariableGet( Mesh % Variables, TRIM(FieldName), ThisOnly )
        ComponentVector = .FALSE.
        
        IF(.NOT. ASSOCIATED(Solution)) THEN
          Solution => VariableGet( Mesh % Variables, TRIM(FieldName)//' 1',ThisOnly )
          IF( ASSOCIATED(Solution)) THEN 
            ComponentVector = .TRUE.
          ELSE
            WRITE(Txt, '(A,A)') 'Nonexistent variable: ',TRIM(FieldName)
            CALL Warn('WriteVtiXMLFile', Txt)
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

