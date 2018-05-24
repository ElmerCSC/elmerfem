!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!/******************************************************************************
! *
! * Subroutine for tracking a particle under the influence of collisions,
! * contacts and external fields.
! *
! ******************************************************************************
! *
! *  Authors: Peter R�back & Juha Ruokolainen
! *  Email:   Peter.Raback@csc.fi & Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 26.5. 2010
! *
! *****************************************************************************/

!> \ingroup{Solvers} 
!> \{

!------------------------------------------------------------------------------
!> The particle-particle and particle-wall interaction for particle dynamics
!> These are placed in a module since a private subroutine may not be passed
!> as a parameter.
!------------------------------------------------------------------------------
MODULE ParticleDynamicsStuff

  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils
  
  IMPLICIT NONE
  
CONTAINS

  !--------------------------------------------------------------------------------------    
  !> Subroutine for getting the force resulting from particle-particle interaction 
  !> This could be used to give forces on granular flow, for example. 
  !--------------------------------------------------------------------------------------    
  SUBROUTINE ParticleParticleContact(dt,Coord,Coord2,Velo,Velo2, &
      Force,Force2, Contact ) 
    
    IMPLICIT NONE

    REAL(KIND=dp):: dt,Coord(3),Coord2(3),Velo(3),Velo2(3),Force(3),Force2(3)
    LOGICAL :: Contact 
    
    REAL(KIND=dp) :: Rad, Mass, Spring, Damping, Friction
    REAL(KIND=dp) :: dist,dr(3),dv(3),eta,rn(3),vn(3),speed,tn1(3),tn2(3)
    REAL(KIND=dp) :: damp_force, spring_force
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found, Visited = .FALSE.
    
    SAVE :: Rad, Mass, Spring, Damping, Friction, Visited

     
    IF( .NOT. Visited ) THEN
      IF( GlobalParticles % NumberOfGroups > 1 ) THEN
        CALL Fatal('ParticleParticleContact','Implemented only for one particle type!')
      END IF
      Params => GetSolverParams()
      Rad = GetCReal(Params,'Particle Radius')
      Mass = GetCReal(Params,'Particle Mass')         
      Spring = GetCReal(Params,'Particle Spring')        
      Damping = GetCReal(Params,'Particle Damping')
      Friction = GetCReal(Params,'Particle Friction',Found)
      Visited =.TRUE.
    END IF
    
    Contact = .FALSE.
    
    ! relative displacement 
    dr = Coord - Coord2
    dist = SQRT( SUM( dr * dr ) ) 
    
    IF( dist < TINY( dist ) ) THEN
      CALL Warn('ParticleParticleContact','Particles are at same point!')
      RETURN
    END IF
    
    ! no contact if the distance is too large
    eta = 2 * Rad - dist
    IF( eta < 0 ) RETURN
    
    ! normal vector     
    rn = dr / dist
    
    ! avoid division by zero at all cost
    dv = Velo - Velo2
    speed = SQRT( SUM( dv * dv ) )
    
    ! if speed is zero, the damping will also be small so no problem with synthetic normal
    IF( speed > TINY( speed ) ) THEN
      vn = dv / speed
    ELSE
      vn = 0.0_dp
      vn(1) = 1.0_dp
    END IF
    
    ! if one needs tangent directions, then activate this
    IF(.FALSE.) THEN
      CALL TangentDirections( rn, tn1, tn2 )
    END IF
    
    ! Currently a linear spring respect to the displacement is given
    ! Here is the place to put the force which could be a complicated function
    ! f=f(r,v,...)
    
    spring_force = eta * Spring 
    !     damp_force = MIN( speed * Damping, spring_force * Friction )
    damp_force = 0.0_dp
    
    
    Force = spring_force * rn - damp_force * vn     
    ! law of force and counterforce:
    Force2 = -Force
    

!     IF( ANY( ISNAN( Force ) ) ) THEN
!       PRINT *,'Force',Force
!       PRINT *,'spring force',spring_force
!       PRINT *,'eta',eta
!       PRINT *,'spring'
!       PRINT *,'vn',vn
!       PRINT *,'damp_force',damp_force
!     END IF

    Contact = .TRUE.
    
  END SUBROUTINE ParticleParticleContact
  

  !---------------------------------------------------------------    
  !> Subroutine for getting particle-particle interaction 
  !> The subroutine may return the new positions and new 
  !> coordinates, or alternatively the initial coordinates are tampered 
  !> so that with standard time-integration the final position will
  !> be the same.
  !---------------------------------------------------------------    
  SUBROUTINE ParticleParticleCollision(dt,Coord,Coord2,Velo,Velo2,&
      Force,Force2, Collision ) 
    
    IMPLICIT NONE

    REAL(KIND=dp):: dt,Coord(3),Coord2(3),Velo(3),Velo2(3),Force(3),Force2(3)
    LOGICAL :: Collision 
    
    REAL(KIND=dp)::  v1na,v2na,v1nb,v2nb
    REAL(KIND=dp) :: Rad1, Rad2, Mass1, Mass2, Coeff 
    REAL(KIND=dp) :: a,b,c,d,dr(3),dv(3),dra(3),rn(3),dta,dtb
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: Found, TrueCollision,SimilarParticles
    LOGICAL :: Visited=.FALSE.
    
    SAVE Visited, SimilarParticles, Rad1, Rad2, Mass1, Mass2, Coeff, &
        TrueCollision
    
    IF(.NOT. Visited ) THEN
      IF( GlobalParticles % NumberOfGroups > 1 ) THEN
        CALL Fatal('ParticleParticleCollision','Implemented only for one particle type')
      END IF
      
      Params => GetSolverParams()
      Rad1 = GetCReal(Params,'Particle Radius',Found)
      IF(.NOT. Found) THEN
        CALL Fatal('ParticleParticleCollision','> Particle Radius < needed!')
      END IF
      Coeff = GetCReal(Params,'Particle Bounciness',Found)
      IF(.NOT. Found ) Coeff = 1.0_dp
      Mass1 = GetCReal(Params,'Particle Mass',Found)  
      Mass2 = Mass1 
      IF(.NOT. Found) THEN
        CALL Fatal('ParticleParticleCollision','> Particle Mass < needed!')
      END IF
      TrueCollision = GetLogical( Params,'True Collision Mode')
      SimilarParticles = .TRUE.
      Visited = .TRUE.
    END IF
    
    Collision = .FALSE.
    
    ! relative displacement and velocity    
    dr = Coord - Coord2
    dv = Velo - Velo2
    
    ! the collision time is found from the conditins |r1(t)-r2(t)|=R1+R2
    ! which results to 2nd order equation for the timestep, here a,b,c
    ! are the coefficicient in the equation. 
    b = SUM( dr * dv )
    
    ! The distance is only growing, there was a collision in history only
    IF( b >= 0.0_dp ) RETURN
    
    a = SUM( dv * dv ) 
    IF( SimilarParticles ) THEN
      c = SUM( dr * dr ) - 4*Rad1**2       
    ELSE
      c = SUM( dr * dr ) - ( Rad1 + Rad2 )**2
    END IF
    d = b*b - a*c
    
    ! negative discriminant means no solution
    IF( d < 0.0_dp ) RETURN
    
    ! time for first collision
    dta = (-b-SQRT(d))/a
    
    ! if time larger than given timestep
    IF( dta >= dt ) RETURN
    
    ! time remaining after the collision
    dtb = dt - dta
    
    ! vector at collision 
    dra = dr + dta * dv
    
    ! normal components at collision 
    rn = dra / SQRT( SUM( dra*dra ) )
    v1na = SUM( Velo * rn )
    v2na = SUM( Velo2 * rn ) 
    
    IF( SimilarParticles ) THEN
      v1nb = ( Coeff * (v2na - v1na) + v1na + v2na ) / 2
      v2nb = ( Coeff * (v1na - v2na) + v2na + v1na ) / 2
    ELSE
      v1nb = ( Coeff * Mass2 * (v2na - v1na) + Mass1 * v1na + Mass2 * v2na ) / ( Mass1 + Mass2 )
      v2nb = ( Coeff * Mass1 * (v1na - v2na) + Mass2 * v2na + Mass1 * v1na ) / ( Mass1 + Mass2 )
    END IF
    
    ! Set either force or velocity directly 
    ! only the normal component of velocity/force is affected by collisions    
    !----------------------------------------------------------------------
    IF( TrueCollision ) THEN
      ! compute the path until the collision
      Coord = Coord + dta * rn * Velo
      Coord2 = Coord2 + dta * rn * Velo2
      
      Velo = Velo + (v1nb-v1na) * rn
      Velo2 = Velo2 + (v2nb-v2na) * rn
      
      ! compute the path after the collision
      Coord = Coord + dtb * Velo
      Coord2 = Coord2 + dtb * Velo2
    ELSE
      Coord = Coord + (v1na-v1nb) * rn * dta
      Coord2 = Coord2 + (v2na-v2nb) * rn * dta
      
      Force = Mass1 * (v1nb-v1na) * rn / dt
      Force2 = Mass2 * (v2nb-v2na) * rn / dt
    END IF
    
    Collision = .TRUE.
    
  END SUBROUTINE ParticleParticleCollision


  !---------------------------------------------------------------    
  !> Subroutine for oding the particle-wall physics at the same time
  !> when locating the particles in the mesh.
  !---------------------------------------------------------------    
  SUBROUTINE ParticleWallProc(No,Rinit,Rfin,Vinit,Vfin,Lambda,FaceIndex,ParticleStatus)

    INTEGER :: No
    REAL(KIND=dp) :: Rinit(3)
    REAL(KIND=dp) :: Rfin(3)
    REAL(KIND=dp) :: Vinit(3)
    REAL(KIND=dp) :: Vfin(3)
    REAL(KIND=dp) :: Lambda
    INTEGER :: FaceIndex
    INTEGER :: ParticleStatus
!------------------------------------------------------------------
    TYPE(Particle_t), POINTER :: Particles
    REAL(KIND=dp) :: UnitVector(3), Normal(3)
    TYPE(Element_t), POINTER :: FaceElement
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: MeshDim
    LOGICAL :: Visited=.FALSE., Found, Stat
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var 
    TYPE(Solver_t), POINTER :: Solver
    REAL(KIND=dp), POINTER :: ParticleEnergy(:), CollisionEnergy(:), &
        ParticleCollisions(:),Basis(:),Reflectivity(:)
    REAL(KIND=dp) :: u,v,w,detJ,TotEnergy,LocalReflectivity,dEnergy
    INTEGER, POINTER :: CollisionPerm(:), NodeIndexes(:),Indexes(:)
    INTEGER :: n,istat,NoParticles
    LOGICAL :: ParticleBunch
    TYPE(ValueList_t), POINTER :: BC

    SAVE Params, Solver, Mesh, MeshDim, ParticleBunch, ElementNodes, Particles, &
        Basis, Reflectivity, Indexes, Visited

    ! when photon reflection visited the first time, do some initialization
    IF ( .NOT. Visited ) THEN 
      Params => GetSolverParams()
      Solver => GetSolver()
      Mesh => GetMesh()
      Particles => GlobalParticles
      MeshDim = Mesh % MeshDim      

      ParticleBunch = GetLogical( Params,'Particle Bunch',Found)

      IF( ParticleBunch ) THEN
        ! get the initial total energy of the photons 
        ! (number of particles times the initial bunch energy)
        ! if it is not found set it 1.0
        !------------------------------------------------------------
        TOTEnergy = ListGetConstReal( Params,'Initial Total Energy',Found)
        IF ( .NOT. found ) TOTEnergy = 1.0_dp
        NoParticles = ListGetInteger( Params,'Number of Particles',Found)
        IF ( .NOT. found ) NoParticles = 1
        
        !------------------------------------------------------------
        CALL ParticleVariableCreate( Particles,'particle energy')
        Var => ParticleVariableGet( Particles,'particle energy')
        ParticleEnergy => Var % Values
        ParticleEnergy = TotEnergy / NoParticles

        !------------------------------------------------------------
        CALL ParticleVariableCreate( Particles,'particle collisions')
        Var => ParticleVariableGet( Particles,'particle collisions')
        ParticleCollisions => Var % Values
        ParticleCollisions = 0

        ! Fetch the Collision energy variable, which is the energy flux 
        ! into the wall provided by the reflected photons
        !---------------------------------------------------------------------
        Var => VariableGet( Mesh % Variables,'Collision Energy')
        IF(.NOT. ASSOCIATED(Var)) THEN
          CALL Info('ParticleDynamics','Creating variable > Collision Energy <')
          CALL VariableAddVector( Mesh % Variables,Mesh,Solver,'Collision Energy')
          Var => VariableGet( Mesh % Variables,'Collision Energy')      
        END IF        
        CollisionEnergy => Var % Values
        CollisionPerm => Var % Perm
      
        n = Mesh % MaxElementNodes
        ALLOCATE( Reflectivity(n), Basis(n), Indexes(n), STAT=istat )
        Reflectivity = 0.0_dp
      
        IF ( istat /= 0 ) THEN
          CALL Fatal( 'ParticleUtils', 'PhotonBunchivity, Memory allocation error' )
        END IF
      END IF
      Visited = .TRUE.
    END IF

    IF( MeshDim == 2 ) THEN
      FaceElement => Mesh % Edges( FaceIndex ) 
    ELSE
      FaceElement => Mesh % Faces( FaceIndex )
    END IF

    ! First advance the particle to the point of collision
    Rinit = Rinit + Lambda * (Rfin - Rinit) 
            
    ! Then reflect the rest assuming fully elastic collision where the 
    ! normal component just switches sign.
    !-----------------------------------------------------------------
    CALL GetElementNodes(ElementNodes, FaceElement )
    Normal = NormalVector( FaceElement, ElementNodes )
    Rfin = Rfin - 2*SUM((Rfin-Rinit)*Normal)*Normal
            
    ! Reorient the velocity vector
    UnitVector = Rfin - Rinit
    UnitVector = UnitVector / SQRT( SUM( UnitVector** 2 ) )
    Vfin = UnitVector * SQRT( SUM( Vinit**2) )    
    
    IF( .FALSE. ) THEN
      PRINT *,'Reflected',No,Lambda
      PRINT *,'Rinit:',Rinit
      PRINT *,'Rfin:',Rfin
      PRINT *,'Vinit:',Vinit
      PRINT *,'Vfin:',Vfin
    END IF

    IF( .NOT. ParticleBunch ) RETURN


    CALL GlobalToLocal( u, v, w, Rinit(1), Rinit(2), Rinit(3), &
        FaceElement, ElementNodes ) 
    stat = ElementInfo( FaceElement, ElementNodes, u, v, w, DetJ, Basis )
    
    ! The BC list is obtained by the fact that that the current element is set    
    BC => GetBC()
    n = FaceElement % TYPE % NumberOfNodes
    NodeIndexes => FaceElement % NodeIndexes
    
    Reflectivity(1:n) = GetReal( BC,'Particle Reflectivity', Stat)
    LocalReflectivity = SUM( Basis(1:n) * Reflectivity(1:n) )

   
    ! Add reflection count and reduce the bunch energy
    !-----------------------------------------------------
    ParticleCollisions(No) = ParticleCollisions(No) + 1
    dEnergy = ( 1 - LocalReflectivity ) * ParticleEnergy(No)
    ParticleEnergy(No) = ParticleEnergy(No) - dEnergy
       
    ! The energy is dumped into the collision node
    !-----------------------------------------------------
    IF( ASSOCIATED( CollisionPerm ) ) THEN
      Indexes(1:n) = CollisionPerm( NodeIndexes(1:n) ) 
    ELSE
      Indexes(1:n) = NodeIndexes(1:n)
    END IF
    CollisionEnergy(Indexes) = CollisionEnergy(Indexes) + dEnergy
  
    
  END SUBROUTINE ParticleWallProc

  
END MODULE ParticleDynamicsStuff


!------------------------------------------------------------------------------
!> Solver particle dynamics equations by utilizing many library functionalities.
!> This solver may take into account interaction between particles. 
!------------------------------------------------------------------------------
SUBROUTINE ParticleDynamics( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils
  USE ParticleDynamicsStuff

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Element_t), POINTER :: CurrentElement
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var
  LOGICAL :: WorkLoop, AssemblyLoop, CollisionInteraction, &
      ContactInteraction, ParticleToField, &
      ParticleInBox, ParticleWall, &
      TrueCollision, StatInfo, ParticleInfo, TimeInfo, Found, &
      FieldReset, ParticlesLocated, ParticleBunch, DoParticleScattering
  INTEGER :: i,j,k,n,dim,NoParticles = 0,&
       ElementIndex, VisitedTimes = 0, nstep, istep, OutputInterval, &
       TimeOrder, TimeStepsTaken=0,estindexes(6),&
       ParticleStepsTaken=0, Group, NoGroups = 0
  REAL(KIND=dp) :: dtime, tottime = 0.0
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: cput1,cput2,dcput
#else
  REAL(KIND=dp) :: cput1,cput2,dcput,CPUTime
#endif
  REAL(KIND=dp), POINTER :: TmpValues(:)
  TYPE(Particle_t), POINTER :: Particles

  SAVE CollisionInteraction, ContactInteraction, NoGroups, &
      ParticleToField, OutputInterval, Nstep, VisitedTimes, &
      TimeOrder, ParticleInBox, &
      tottime, TimeStepsTaken, ParticleStepsTaken,ParticleWall,StatInfo, &
      ParticleInfo, TimeInfo, TrueCollision, FieldReset, &
      ParticlesLocated,DoParticleScattering


!------------------------------------------------------------------------------

  CALL Info('ParticleDynamics','-----------------------------------------', Level=4 )
  CALL Info('ParticleDynamics','Following the path of the particles',Level=4) 

  VisitedTimes = VisitedTimes + 1

  Particles => GlobalParticles
  PSolver => Solver
  Params => GetSolverParams()
  Mesh => Solver % Mesh
  DIM = CoordinateSystemDimension()

  ! Do some initalialization: allocate space, check fields  
  !------------------------------------------------------------------------
  IF( VisitedTimes == 1 ) THEN
    TimeOrder = GetInteger( Params,'Time Order',Found)
    IF(.NOT. Found) TimeOrder = 2

    NoGroups = GetInteger( Params,'Number Of Particle Groups',Found )
    IF( Found ) THEN
      Particles % NumberOfGroups = NoGroups 
    ELSE
      ! This means that group concept is passive
      ! We want one group to be already a test case for the group concept
      NoGroups = 0
    END IF
      
    CALL SetParticlePreliminaries( Particles, dim, TimeOrder )

    i = GetInteger( Params,'Random Seed',Found ) 
    IF( Found ) CALL RANDOM_SEED(i)

    ParticleToField = GetLogical( Params,'Particle To Field',Found)
    CollisionInteraction = GetLogical( Params,'Particle Particle Collision',Found)
    ContactInteraction = GetLogical( Params,'Particle Particle Contact',Found)
    OutputInterval = GetInteger( Params,'Output Interval',Found)

    Nstep = GetInteger( Params,'Max Timestep Intervals',Found)
    IF(.NOT. Found) Nstep = 1
    ParticleInBox = GetLogical( Params,'Box Particle Collision',Found) .OR. &
        GetLogical( Params,'Box Particle Contact',Found)

    ParticleWall = GetLogical( Params,'Wall Particle Contact',Found) 

    DoParticleScattering = GetLogical( Params,'Particle Scattering',Found) 

    TrueCollision = .FALSE.
    CALL ListAddLogical( Params,'True Collision Mode',TrueCollision)

    StatInfo = GetLogical( Params,'Statistical Info',Found)
    ParticleInfo = GetLogical( Params,'Particle Info',Found)
    TimeInfo = GetLogical( Params,'Timing Info',Found)
    ParticlesLocated = .FALSE.
    
    IF( ParticleToField ) THEN
      FieldReset = GetLogical( Params,'Particle To Field Reset',Found)
    END IF
  END IF


  ! Initialize particles at first time visiting, or each time if requested
  !-------------------------------------------------------------------------
  IF( VisitedTimes == 1 .OR. &
      GetLogical( Params,'Reinitialize Particles',Found) ) THEN

    IF( NoGroups > 0 ) THEN
      IF( ListGetLogical( Params,'Set Particle Group By Domain',Found ) ) THEN
        CALL InitializeParticles( Particles, Group = 1 )
        CALL Info('ParticleDynamics','Setting particle group by domain',Level=9)
        DO i=1, Particles % NumberOfParticles
          j = Particles % ElementIndex(i)
          IF( j == 0 ) THEN
            CALL LocateParticles(Particles, ParticleWallProc ) 
            j = Particles % ElementIndex(i)
          END IF
          Particles % Group(i) = Mesh % Elements(j) % BodyId 
        END DO
        DO Group = 1, NoGroups
          j = COUNT( Particles % Group == Group )
          CALL Info('ParticleDynamics','Group '//TRIM(I2S(Group))//' particles: '//TRIM(I2S(j)),Level=9)
        END DO
      ELSE
        DO Group = 1, NoGroups
          CALL Info('ParticleDynamics','Initializing particles in group '//TRIM(I2S(group)),Level=5)
          CALL ListPushNameSpace('group'//TRIM(I2S(Group))//':')
          CALL InitializeParticles( Particles, AppendParticles = .TRUE.,Group = Group )
          CALL ListPopNameSpace()
        END DO
      END IF
    ELSE
      CALL InitializeParticles( Particles )
    END IF
    
    ParticlesLocated = .FALSE.
    IF( GetLogical( Params,'Particle Distance',Found) ) THEN
      CALL ParticleVariableCreate( Particles,'Particle Distance' )
    END IF
  END IF

  CALL ReleaseWaitingParticles(Particles)   

  IF(.FALSE. .AND. StatInfo) THEN
    CALL ParticleStatistics( Particles, 0 )
    CALL ParticleStatistics( Particles, 1 )
  END IF

  ! a logaritmic scale of indexes is used to estimate time
  !--------------------------------------------------------
  IF( TimeInfo ) THEN
    cput1 = CPUTime()
    estindexes(1) = nstep / 2
    DO i=1,5
      estindexes(i+1) = estindexes(i) / 10
    END DO
  END IF


  DO i=1,nstep

    ! Get the timestep size, initialize at 1st round
    !--------------------------------------------------------------
    dtime = GetParticleTimeStep( Particles, i == 1 )
!.NOT. ParticlesLocated )

    IF( ParticleToField .AND. i > 1) THEN      
      AssemblyLoop = .NOT. FieldReset
    ELSE
      AssemblyLoop = .FALSE.
    END IF
    
    ! If size of timestep goes to zero then no more steps are needed
    !---------------------------------------------------------------
    IF( dtime < TINY( dtime ) ) EXIT
    tottime = tottime + dtime


    ! Initialize the timestep, in practive just set force to zero
    !--------------------------------------------------------------
    CALL ParticleInitializeTime( Particles )

    ! If there are periodic BCs apply them just before locating the particles
    !------------------------------------------------------------------------
    ! CALL ParticleBoxPeriodic( Particles )

    ! Find the elements (and only the elements) in which the particles are in 
    ! This might also not be done for some particle-only problems but those 
    ! are probably not relevant in conjunction with Elmer.
    !------------------------------------------------------------------------
    IF( .NOT. ( i==1 .AND. ParticlesLocated ) ) THEN
      CALL LocateParticles(Particles, ParticleWallProc ) 
      ParticlesLocated = .TRUE.
    END IF

    ! Eliminate the particles that sit on the wall but are not on
    ! partition interface.
    !------------------------------------------------------------------------
    CALL EliminateExitingParticles( Particles )

    ! Calculate the force resulting from external fields in mesh
    ! and / or cumpulate the r.h.s. of matrix equation with data 
    !------------------------------------------------------------------
    CALL ParticleFieldInteraction( Particles, dtime, .TRUE. , AssemblyLoop ) 
    
    NoParticles = Particles % NumberOfParticles
    ParticleStepsTaken = ParticleStepsTaken + NoParticles
    TimeStepsTaken = TimeStepsTaken + 1
    
    ! Interaction with the walls
    !---------------------------------------------------------------
    IF( ParticleInBox ) THEN
      ! Faster for rectangular and hexahedral domains
      !---------------------------------------------------------------
      CALL  ParticleBoxContact( Particles ) 
    ELSE IF( ParticleWall ) THEN
      ! Generic version
      !---------------------------------------------------------------
      CALL ParticleWallContact( Particles, dtime )
    END IF
    
    ! If there is either collisions or contacts between particles
    ! create the structures for closest neighbours
    !------------------------------------------------------------
    IF( CollisionInteraction .OR. ContactInteraction ) THEN
      CALL CreateNeighbourList( Particles )

      ! Add the forces in case of collision contacts where the
      ! effect is mutated into acceleration values.
      !------------------------------------------------------------
      IF( CollisionInteraction ) THEN
        CALL ParticleParticleInteraction( Particles, dtime, .TRUE., &
            ParticleParticleCollision ) 
      END IF

      ! Add contact interaction, for example due to granular forces
      !------------------------------------------------------------
      IF( ContactInteraction ) THEN
        CALL ParticleParticleInteraction( Particles, dtime, .FALSE., &
            ParticleParticleContact ) 
      END IF
      
      ! In parallel case destroy the ghost particles needed for 
      !---------------------------------------------------------------
      CALL DestroyGhostParticles( Particles ) 
    END IF

    ! Do the update for particle velocities and positions
    ! v = v0 + at, r = r0 + vt
    !---------------------------------------------------------------
    CALL ParticleAdvanceTimestep( Particles )
    
    ! Do particle scattering from the bulk, for example acoustic scattering
    !----------------------------------------------------------------------
    IF( DoParticleScattering ) THEN
      CALL ParticleScattering( Particles ) 
    END IF
    
    ! If there are periodic BCs apply them just before locating the particles
    !------------------------------------------------------------------------
    CALL ParticleBoxPeriodic( Particles )

    ! Delete the particles that were not found anymore
    !---------------------------------------------------------------
    CALL DeleteLostParticles( Particles )

    IF( OutputInterval > 0 .AND. MOD(i,OutputInterval) == 0) THEN
      CALL SaveParticleData( Model,Solver,dt,TransientSimulation )
    END IF

    ! Write estimates of remaining time in log scale
    !---------------------------------------------------------------
    IF( TimeInfo) THEN
      IF( ANY( estindexes == i ) ) THEN
        cput2 = CPUTime()
        dcput = cput2 - cput1
        IF( dcput > 0.5 ) THEN
          WRITE( Message,'(A,F8.3)') 'Fraction computed (s)  :',(100.0_dp)*i/nstep
          CALL Info('ParticleDynamics',Message)
          WRITE( Message,'(A,F8.3)') 'Consumed time (s)  :',dcput
          CALL Info('ParticleDynamics',Message)
          WRITE( Message,'(A,F8.3)') 'Remaining time (s) :',dcput*(nstep-i)/i
          CALL Info('ParticleDynamics',Message)
        END IF
      END IF
    END IF

  END DO

  ! Do one last assembly if particle field is requested
  !------------------------------------------------------------------------
  IF( ParticleToField ) THEN
!    CALL ParticleBoxPeriodic( Particles )
    CALL LocateParticles(Particles, ParticleWallProc ) 
    CALL EliminateExitingParticles( Particles )
    CALL ParticleFieldInteraction( Particles, dtime, .FALSE., .TRUE. )
    ParticlesLocated = .TRUE.
  END IF



  ! In the end, compute the fields
  ! Interaction with the fields is typically with external solvers
  ! so no idea to do it after each timestep.
  !---------------------------------------------------------------   
  IF(StatInfo) THEN
    CALL ParticleStatistics( Particles, 0 )
    CALL ParticleStatistics( Particles, 1 )
  END IF
  
  IF( ParticleInfo ) THEN
    CALL ParticleInformation(Particles, ParticleStepsTaken, &
	TimeStepsTaken, tottime )
  END IF
    
  CALL Info('ParticleDynamics','All done',Level=4)
  CALL Info('ParticleDynamics', '-----------------------------------------', Level=4 )
  
  
CONTAINS   
  
  !---------------------------------------------------------
  !> Advance the particles with a time step. The timestep may
  !> also be an intermediate Runge-Kutta step.
  !---------------------------------------------------------
  SUBROUTINE ParticleScattering( Particles )
    TYPE(Particle_t), POINTER :: Particles
    
    REAL(KIND=dp) :: dt,dt0,Angle,VeloAbs,Pscatter,&
        Velo(3),Ptest,MFP
    INTEGER :: No,NoParticles,NoScattered,Status
    LOGICAL :: Visited=.FALSE.,Scatter,Found
    
    SAVE Visited,MFP
    
    IF(.NOT. Visited ) THEN
      Params => GetSolverParams()
      TimeOrder = Particles % TimeOrder
      dim = Particles % dim
      
      IF( .NOT. Particles % DtConstant ) THEN
        CALL Fatal('ParticleScattering','Timestep should not be variable in this routine!')
      END IF
      
      MFP = ListGetConstReal( Params,'Particle MFP',Found ) 
      IF(.NOT. Found ) THEN
        CALL Fatal('ParticleScattering','Keyword > Particle MFP < not given!')
      END IF
      
      Visited = .TRUE.
    END IF
    
    dt = Particles % dTime
    NoParticles = Particles % NumberOfParticles
    NoScattered = 0
    
    DO No=1, NoParticles
      
      Status = Particles % Status(No)
      
      IF ( Status >= PARTICLE_LOST ) CYCLE
      IF ( Status <= PARTICLE_INITIATED ) CYCLE
      IF ( Status == PARTICLE_WALLBOUNDARY ) CYCLE
      
      Velo(1:dim) = Particles % Velocity(No,:) 	
      VeloAbs = SQRT( SUM( Velo(1:dim) ** 2 ) )
      
      Pscatter = 1.0_dp - EXP( -(Dt*VeloAbs) / MFP ) 
      Ptest = EvenRandom()
      
      Scatter = ( Ptest < Pscatter )  
      
      IF( Scatter ) THEN        
        Angle = 2*PI*EvenRandom()        
        Particles % Velocity(No,1) = COS( Angle ) * Velo(1) - SIN( Angle ) * Velo(2)
        Particles % Velocity(No,2) = COS( Angle ) * Velo(2) + SIN( Angle ) * Velo(1)

        ! For higher order scheme enforce the previous velocity also to zero since
        ! otheriwse there is a funny velocity correction added.
        IF( Particles % TimeOrder > 1 ) THEN
          Particles % PrevVelocity(No,:) = Particles % Velocity(No,:)
        END IF
        NoScattered = NoScattered + 1
      END IF
    END DO
    
    WRITE(Message,'(A,I0)') 'Number of Scattered particles: ',NoScattered
    CALL Info('ParticleScattering',Message,Level=9)
    
    CALL ListAddConstReal( Model % Simulation,'res: Scattered particles',&
        1.0_dp*NoScattered)

  END SUBROUTINE ParticleScattering
  




  !------------------------------------------------------------------------
  ! Compute field values at the given points in the FE mesh. 
   !-------------------------------------------------------------------------
   SUBROUTINE ParticleFieldInteraction(Particles,dtime,SetParticles,SetFields ) 
     
     TYPE(Particle_t), POINTER :: Particles
     REAL(KIND=dp) :: dtime 
     LOGICAL :: SetParticles, SetFields
     !-------------------------------------------------------------------------

     TYPE VarPointer_t
       TYPE(Variable_t), POINTER :: Var
     END TYPE VarPointer_t

#define MAXPARFIELDS 20
     TYPE(VarPointer_t) :: ActiveVars(MAXPARFIELDS)
     TYPE(Element_t), POINTER :: BulkElement
     INTEGER :: No, Status
     REAL(KIND=dp) :: Coord(3),Velo(3), Force(3)
     
     TYPE(Element_t), POINTER :: BulkElement2
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(Valuelist_t), POINTER :: Params
     REAL(KIND=dp) :: PotAtPoint, GradPotAtPoint(3),VeloAtPoint(3), &
         GradVeloAtPoint(3,3)
     LOGICAL :: Stat, UseGradVelo, CoordCond, VeloCond, Visited = .FALSE., &
         GotIt, GotPot, GotPot2, GotVelo
     INTEGER :: i,j,k,l,n,dim,TimeOrder,NoGroups, MaxField, PrevGroup, CurrGroup
     INTEGER, POINTER :: NodeIndexes(:)
     REAL(KIND=dp) :: SqrtElementMetric, Weight, TimeDecay, DistDecay, Dist, &
         ParticleVolume, FluidDensity, VolumeFraction,UserCoeff,&
         Cchar,Tchar,Prevdtime
     REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
     REAL(KIND=dp) :: val, ValCoeff, Gravity(3), sumf, sumw
     REAL(KIND=dp), POINTER :: gWork(:,:), ForceVector(:)
     INTEGER, POINTER :: ForcePerm(:)
     CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, DensityName, FieldMode, FieldWeight, GroupName, str
     TYPE(Variable_t), POINTER :: VeloVar, PotVar, PotVar2, VeloCondVar, CoordCondVar, WeightVar
     LOGICAL :: GotGravity, GotBuoyancy, GotField, &
         GotTimeDecay, GotDistDecay, GotFieldMode, GotFieldWeight, &
         NormalizedVars(MAXPARFIELDS)
     INTEGER :: ActiveOpers(MAXPARFIELDS),ActiveGroups(MAXPARFIELDS)
     TYPE(Variable_t), POINTER :: DistVar
     
     REAL(KIND=dp) :: mass, damping, charge, dragcoeff, rad
     REAL(KIND=dp), POINTER :: massv(:), dampingv(:), chargev(:), dragcoeffv(:), radv(:)
     LOGICAL :: GotMass, GotDamping, GotCharge, GotDrag, GotRad
     
     SAVE :: Visited, dim, Basis, dBasisdx, &
         FieldMode, FieldWeight, TimeDecay, DistDecay, UseGradVelo, TimeOrder, &
         GotFieldMode, GotFieldWeight, GotGravity, GotDamping, GotTimeDecay, GotDistDecay, &
         GotPot, GotPot2, GotVelo, Gravity, Damping, VeloCond, CoordCond, GotBuoyancy, &
         ParticleVolume, GotField, CoordCondVar, VeloCondVar, DensityName, &
         PotVar, PotVar2, VeloVar, Mesh, PrevDtime, DistVar, &
         ActiveVars, ActiveOpers, ActiveGroups, NormalizedVars, MaxField, NoGroups



     Params => GetSolverParams()

     IF( .NOT. Visited ) THEN
       Mesh => GetMesh()
       dim = Mesh % MeshDim
       n = Mesh % MaxElementNodes
       ALLOCATE( Basis(n), dBasisdx(n, 3) )

         
       GotBuoyancy = GetLogical( Params,'Particle Lift',Found)

       GotGravity = GotBuoyancy .OR. ListGetLogical( Params,'Particle Gravity',Found)
       IF( GotGravity ) THEN
         gWork => ListGetConstRealArray( CurrentModel % Constants,'Gravity',Found)
         IF ( ASSOCIATED(gwork) ) THEN
           Gravity = gWork(4,1) * gWork(1:3,1)
         ELSE
           CALL Fatal('ParticleFieldInteraction','Gravity and Lift requires gravity!')
         END IF
       END IF

       VariableName = ListGetString(Params,'Potential Variable Name',GotPot)
       IF( GotPot ) THEN
         PotVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
         IF(.NOT. ASSOCIATED( PotVar ) ) THEN
           CALL Fatal('ParticleFieldInteraction','Potential field variable does not exist: '//TRIM(VariableName))           
         END IF
       END IF
       
       VariableName = ListGetString(Params,'Secondary Potential Variable Name',GotPot2)
       IF( GotPot2 ) THEN
         PotVar2 => VariableGet( Mesh % Variables, TRIM(VariableName) )
         IF(.NOT. ASSOCIATED( PotVar2 ) ) THEN
           CALL Fatal('ParticleFieldInteraction','Potential field variable does not exist: '//TRIM(VariableName))           
         END IF
       END IF

       VariableName = ListGetString(Params,'Velocity Variable Name',GotVelo)
       IF( GotVelo ) THEN
         VeloVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
         IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
           CALL Fatal('ParticleFieldInteraction','Velocity field variable does not exist: '//TRIM(VariableName))           
         END IF
         UseGradVelo = GetLogical( Params,'Velocity Gradient Correction',Found)
       ELSE
         UseGradVelo = .FALSE.
       END IF
    
       VariableName = ListGetString(Params,'Velocity Condition Variable Name',VeloCond)
       IF( VeloCond ) THEN
         VeloCondVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
         IF(.NOT. ASSOCIATED( VeloCondVar ) ) THEN
           CALL Fatal('ParticleFieldInteraction','Velocity condition field variable does not exist: '//TRIM(VariableName))           
         END IF                  
       END IF

       VariableName = ListGetString(Params,'Coordinate Condition Variable Name',CoordCond)
       IF( CoordCond ) THEN
         CoordCondVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
         IF(.NOT. ASSOCIATED( CoordCondVar ) ) THEN
           CALL Fatal('ParticleFieldInteraction','Coordinate condition field variable does not exist: '//TRIM(VariableName))           
         END IF                  
       END IF

       TimeOrder = Particles % TimeOrder

       IF( ParticleToField ) THEN
         FieldWeight = GetString( Params,'Particle To Field Weight',GotFieldWeight)
         DistDecay = GetCReal( Params,'Particle Decay Distance',GotDistDecay)
         TimeDecay = GetCReal( Params,'Particle Decay Time',GotTimeDecay)         

         FieldMode = GetString( Params,'Particle To Field Mode',GotFieldMode)
         IF( GotFieldMode ) THEN
           CALL Fatal('ParticleFieldInteraction','> Particle to Field Mode < is obsolite. Use > Field 1 < instaed!')
         END IF
           
         ActiveOpers = 0
         ActiveGroups = 0
         NormalizedVars = .FALSE.

         i = 0
         DO i=1,MAXPARFIELDS

           WRITE( str,'(A,I0)') 'Field ',i
           FieldMode = ListGetString( Params, str, GotIt )
           IF( .NOT. GotIt ) THEN
             MaxField = i - 1             
             EXIT          
           END IF

           SELECT CASE( FieldMode ) 

           CASE ('weight')             
             j = 1
           CASE ('energy') 
             j = 2
           CASE ('kinetic energy')
             j = 3
           CASE ('potential energy')
             j = 4
           CASE ('electrostatic energy')
             j = 5
           CASE('charge')
             j = 6
           CASE ('speed') 
             j = 7
           CASE ('force') 
             j = 8
           CASE DEFAULT
             CALL Fatal('ParticleFieldInteraction','Unknown field mode: '//TRIM(FieldMode))

           END SELECT

           ActiveOpers(i) = j

           j = 0
           IF( Particles % NumberOfGroups > 0 ) THEN
             WRITE( str,'(A,I0)') 'Group ',i
             j = ListGetInteger( Params, str, GotIt )             
             ActiveGroups(i) = j
           END IF

           IF( j > 0 ) THEN
             WRITE( str,'(A,I0)') 'Group Name ',i
             GroupName = ListGetString( Params, str, GotIt)
             IF( .NOT. GotIt ) CALL Fatal('ParticleDynamics','Unfound keyword: '//TRIM(str))
             VariableName = 'Particle '//TRIM(FieldMode)//' '//TRIM(GroupName)
           ELSE
             VariableName = 'Particle '//TRIM(FieldMode)
           END IF
             
           Var => VariableGet( Mesh % Variables,VariableName )
           IF(.NOT. ASSOCIATED(Var)) THEN
             CALL Info('ParticleDynamics','Creating variable: '//VariableName ) 
             CALL VariableAddVector( Mesh % Variables,Mesh,PSolver,VariableName )
             Var => VariableGet( Mesh % Variables,VariableName)      
           END IF

           ActiveVars(i) % Var => Var

           WRITE( str,'(A,I0)') 'Field Normalize ',i         
           NormalizedVars(i) = ListGetLogical( Params, str, GotIt )
         END DO

         CALL Info('ParticleDynamics','Number of particle fields: '//TRIM(I2S(MaxField)),Level=6)
       END IF
       
       DensityName = 'Density'

       IF( GotDistDecay ) THEN
         DistVar => ParticleVariableGet( Particles,'Particle Distance' )
         IF( .NOT. ASSOCIATED( DistVar ) ) THEN
           CALL Fatal('ParticleDynamics','> Particle Distance < should exist as variable!')
         END IF
       END IF

       GotField = ParticleToField .OR. GotVelo .OR. GotPot .OR. &
           VeloCond .OR. CoordCond .OR. GotBuoyancy

       Visited = .TRUE.
     END IF

     ! Normalization and resetting of particle effect on fields
     !-------------------------------------------------------------------------     
     IF( ParticleToField ) THEN       
       UserCoeff = GetCReal( Params,'Particle To Field Coefficient',Found)
       IF(.NOT. Found) UserCoeff = 1.0_dp

       Tchar = GetCReal( Params,'Field Decay Time',Found)
       IF( Found ) THEN
         Cchar = EXP( -dtime / Tchar ) 
         DO i=1,MaxField
           ActiveVars(i) % Var % Values = Cchar * ActiveVars(i) % Var % Values
         END DO
       END IF
              
       IF( FieldReset ) THEN
         DO i=1,MaxField
           ActiveVars(i) % Var % Values = 0.0_dp
         END DO
       END IF
     END IF

     NoParticles = Particles % NumberOfParticles
     NoGroups = Particles % NumberOfGroups     
     PrevGroup = -1


     ! The many groups case is treated separately since it adds limitation to the keywords being
     ! constant. For one group the parameters could depend on global parameters such as time. 
     !------------------------------------------------------------------------------------------
     IF( NoGroups > 1 ) THEN
       massv => ListGetConstRealArray1( Params,'Particle Mass',GotMass)
       mass = 0.0_dp
       dampingv => ListGetConstRealArray1( Params,'Particle Damping',GotDamping)
       damping = 0.0_dp
       Radv => ListGetConstRealArray1(Params,'Particle Radius',GotRad )
       rad = 0.0_dp
       chargev => ListGetConstRealArray1( Params,'Particle Charge',GotCharge)
       charge = 0.0_dp
       dragcoeffv => ListGetConstRealArray1( Params,'Particle Drag Coefficient',GotDrag)
       dragcoeff = 0.0_dp
     ELSE
       mass = GetCReal( Params,'Particle Mass',GotMass)
       damping = GetCReal( Params,'Particle Damping',GotDamping)         
       Rad = GetCReal(Params,'Particle Radius',GotRad)
       charge = GetCReal( Params,'Particle Charge',GotCharge)
       dragcoeff = GetCReal( Params,'Particle Drag Coefficient',GotDrag)
     END IF
         
     IF( GotBuoyancy .AND. .NOT. GotRad ) THEN
       CALL Fatal('ParticleFieldInteraction','> Particle Radius < is needed for buoyancy!')
     END IF
     
     IF( GotGravity .AND. .NOT. GotMass ) THEN
       CALL Warn('ParticleFieldInteraction','> Particle Mass < is needed by gravity!')
     END IF

     IF( GotPot .AND. .NOT. GotCharge) THEN
       CALL Fatal('ParticleFieldInteraction',&
           '> Particle Charge < is needed by external field!')
     END IF
      
     IF( GotVelo .AND. .NOT. GotDrag ) THEN
       CALL Fatal('ParticleFieldInteraction','> Particle Drag Coefficient < required with velocity!')        
     END IF

     IF( GotBuoyancy ) THEN
       IF( dim == 2 ) THEN
         ParticleVolume = PI * Rad ** 2
       ELSE
         ParticleVolume = (4.0_dp/3) * PI * Rad ** 3
       END IF
     END IF
       
     
     
     DO No = 1, NoParticles

       Status = GetParticleStatus( Particles, No )

       
       IF( Status >= PARTICLE_LOST ) CYCLE
       IF( Status <= PARTICLE_INITIATED ) CYCLE

       
       ElementIndex = GetParticleElement( Particles, No )

       IF( ElementIndex == 0 ) CYCLE

       BulkElement => Mesh % Elements( ElementIndex )
       IF(.NOT. ASSOCIATED( BulkElement ) ) CYCLE

       Coord = GetParticleCoord( Particles, No )
       Velo = GetParticleVelo( Particles, No )
       Force = 0.0_dp       
       
       IF( NoGroups > 1 ) THEN
         CurrGroup = GetParticleGroup( Particles, No )        
         
         IF( CurrGroup /= PrevGroup ) THEN
           IF(GotMass) mass = massv(MIN(SIZE(massv),CurrGroup))
           IF(GotDamping) damping = dampingv(MIN(SIZE(dampingv),CurrGroup))
           IF(GotRad) rad = radv(MIN(SIZE(radv),CurrGroup))
           IF(GotCharge) charge = chargev(MIN(SIZE(chargev),CurrGroup))
           IF(GotDrag) dragcoeff = dragcoeffv(MIN(SIZE(dragcoeffv),CurrGroup))

           IF( GotBuoyancy ) THEN
             IF( dim == 2 ) THEN
               ParticleVolume = PI * Rad ** 2
             ELSE
               ParticleVolume = (4.0_dp/3) * PI * Rad ** 3
             END IF
           END IF
           
           PrevGroup = CurrGroup
         END IF
       END IF
              
       !-------------------------------------------------------------------------
       ! Add constant fields i.e. gravity and constant damping force
       !-------------------------------------------------------------------------
       IF( GotGravity ) THEN
         Force = Force + Gravity * Mass
       END IF

       ! For 1st order models the velocity is solved implicitely, when drag is known
       ! Therefore it is a explicit force only for 2nd order models.
       !----------------------------------------------------------------------------
       IF( GotDamping .AND. TimeOrder == 2 ) THEN
         Force = Force - damping * Velo * Mass 
       END IF
              
       IF( GotField ) THEN       
         IF( GotPot .OR. UseGradVelo ) THEN
           stat = ParticleElementInfo( BulkElement, Coord, &
               SqrtElementMetric, Basis, dBasisdx )
         ELSE
           stat = ParticleElementInfo( BulkElement, Coord, &
               SqrtElementMetric, Basis )
         END IF
         IF( .NOT. stat ) THEN
           CALL Warn('ParticleFieldInteraction','Particle not in element')
           CYCLE
         END IF
         
         !-------------------------------------------------------------------------
         ! Set Dirichlet conditions for velocity / coordinate
         ! The condition is computed from an external field and hence this is the 
         ! appropriate place to set this flag.
         !-------------------------------------------------------------------------
         IF( VeloCond ) THEN
           CALL GetScalarFieldInMesh(VeloCondVar, BulkElement, Basis, val ) 
           IF( val > TINY( val ) ) Status = PARTICLE_FIXEDVELO
         END IF

         IF( CoordCond ) THEN
           CALL GetScalarFieldInMesh(CoordCondVar, BulkElement, Basis, val ) 
           IF( val > TINY( val ) ) Status = PARTICLE_FIXEDCOORD
         END IF
         
         !-------------------------------------------------------------------------
         ! Add interaction with fields i.e. fluidic, electrostatic forces etc.
         ! If acceleration is not needed solve implicitely for the velocity
         ! when drag coefficient is known. Hence do not set that as force here.
         !-------------------------------------------------------------------------
         IF( GotVelo ) THEN
           IF( UseGradVelo ) THEN       
             CALL GetVectorFieldInMesh(VeloVar,BulkElement, Basis, VeloAtPoint, &
                 dBasisdx, GradVeloAtPoint )
             DO i=1,dim
               VeloAtPoint(i) = VeloAtPoint(i) + &
                   0.5_dp * SUM( GradVeloAtPoint(i,1:dim) * Velo(1:dim) ) * dtime        
             END DO
           ELSE
             CALL GetVectorFieldInMesh(VeloVar, BulkElement, Basis, VeloAtPoint )
           END IF
           
           IF( TimeOrder == 2 ) THEN
             Force = Force + dragcoeff * ( VeloAtPoint - Velo )  
           ELSE     
             Force = Force + dragcoeff * VeloAtPoint 
           END IF
         ELSE
           VeloAtPoint = 0.0_dp
         END IF

         IF( GotPot ) THEN
           CALL GetScalarFieldInMesh(PotVar, BulkElement, Basis, PotAtPoint, &
               dBasisdx, GradPotAtPoint )
           Force = Force - charge * GradPotAtPoint 
         END IF
         
         ! there can be a secondary potential field also
         IF( GotPot2 ) THEN
           CALL GetScalarFieldInMesh(PotVar2, BulkElement, Basis, PotAtPoint, &
               dBasisdx, GradPotAtPoint )
           Force = Force - charge * GradPotAtPoint 
         END IF

         IF( GotBuoyancy ) THEN
           IF( GetParticleElementIntersection( Particles, BulkElement, Basis, Coord, &
               Rad, BulkElement2, VolumeFraction ) ) THEN
             FluidDensity = GetMaterialPropertyInMesh(DensityName, BulkElement, Basis, &
                 BulkElement2, VolumeFraction ) 
           ELSE
             FluidDensity = GetMaterialPropertyInMesh(DensityName, BulkElement, Basis )
           END IF
           Force = Force - Gravity * ParticleVolume * FluidDensity
         END IF
         
         !-------------------------------------------------------------------------
         ! the value at point is obtained from a property of the particles
         ! which may be accumulated with time. Note that the weight could
         ! be also ~1/r^2 from the nodes etc. 
         !-------------------------------------------------------------------------
       END IF

       
       IF( SetFields ) THEN
         NodeIndexes =>  BulkElement % NodeIndexes
         n = BulkElement % Type % NumberOfNodes

         ! Weight depending on the particle tracking resolution
         !------------------------------------------------------
         IF( GotFieldWeight ) THEN
           
           SELECT CASE( FieldWeight ) 
             
           CASE ( 'distance' )
             ValCoeff = SQRT( SUM( Velo ** 2 ) ) * prevdtime
             
           CASE ( 'time' ) 
             ValCoeff = prevdtime
             
           CASE( 'speed' ) 
             ValCoeff = SQRT( SUM( Velo ** 2 ) )
             
           CASE DEFAULT 
             CALL Fatal('ParticleFieldInteraction','Unknown field weight: '//TRIM(FieldMode))
             
           END SELECT
         ELSE
           ValCoeff = 1.0_dp
         END IF

         ! Decay of particles weight in time or space
         !---------------------------------------------------
         IF( GotTimeDecay ) THEN
           ValCoeff = ValCoeff * EXP(-tottime/TimeDecay)            
         END IF
         
         IF( GotDistDecay ) THEN
           dist = DistVar % Values( No ) 
           ValCoeff = ValCoeff * EXP(-dist/DistDecay)
         END IF
         
         ! User defined normalization
         !---------------------------------------------------
         ValCoeff = UserCoeff * ValCoeff

         DO k=1,MaxField
           l = ActiveGroups(k)
           IF( l > 0 ) THEN
             IF( CurrGroup /= l ) CYCLE
           END IF
          
           ! physical weight dependent of the particle configuration
           !--------------------------------------------------------
           SELECT CASE( ActiveOpers(k) ) 

           CASE( 1 ) 
             val = 1.0_dp

           CASE( 2 )
             ! ( kinetic + potential + electrostatic ) energy
             val = 0.5 * Mass * SUM( Velo(1:dim)**2 ) + &
                 Mass * SUM( Gravity(1:dim) * Coord(1:dim) ) + &
                 Charge * PotAtPoint

           CASE( 3 )
             val = 0.5 * Mass * SUM( Velo ** 2 ) 

           CASE( 4 ) 
             val = Mass * SUM( Gravity(1:dim) * Coord(1:dim) )

           CASE( 5 )
             val = Charge * PotAtPoint

           CASE( 6 ) 
             val = Charge 

           CASE( 7 )
             val = SQRT( SUM( VeloAtPoint ** 2 ) )

           CASE( 8 ) 
             val = SQRT( SUM( Force(1:dim) ** 2 ) )

           END SELECT

           val = val * ValCoeff 
           Var => ActiveVars(k) % Var
           ForceVector => Var % Values

           IF( ASSOCIATED( Var % Perm ) ) THEN
             ForcePerm => Var % Perm 
           ELSE
             ForcePerm => NULL()
           END IF                          

           DO i = 1,n
             ! As the weight should be proporpotional to the particle amount rather than
             ! element volume this is not to be multiplied with local element size!
             !--------------------------------------------------------------------------
             weight = Basis(i)
             
             j = NodeIndexes(i) 
             IF( ASSOCIATED( ForcePerm ) ) THEN
               j = ForcePerm( j )
               IF( j == 0 ) CYCLE
             END IF
             
             ForceVector( j ) = ForceVector( j ) + weight * val
           END DO

         END DO
       END IF

       IF( SetParticles ) THEN
         CALL AddParticleForce( Particles, No, Force )
         CALL SetParticleStatus( Particles, No, Status )             
       END IF
       
     END DO

     IF( SetFields ) THEN
       VariableName = ListGetString(Params,'Particle Nodal Weights',GotIt)
       IF( GotIt ) THEN
         CALL Info('ParticleFieldInteraction','Setting the average to zero',Level=9)
         WeightVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
         sumw = SUM( WeightVar % Values )
         
         DO j=1,MAXPARFIELDS
           IF( NormalizedVars(j) ) THEN
             ForceVector => ActiveVars(j) % Var % Values
             sumf = SUM( ForceVector )
             IF( SIZE( WeightVar % Values ) /= SIZE( ForceVector ) ) THEN
               CALL Fatal('ParticleFieldInteraction','Sizes are assumed to be the same')
             END IF
             
             IF( ABS( sumf ) > TINY( sumf ) ) THEN       
               sumw = SUM( WeightVar % Values )
               ForceVector = ForceVector - (sumf / sumw) * WeightVar % Values
             END IF
           END IF
         END DO
       END IF
     END IF

     PrevDtime = dtime
     

   END SUBROUTINE ParticleFieldInteraction
   

 
   
   !---------------------------------------------------------------    
   ! Checks the boundaries for general limits. 
   ! The radius is still assumed to be constant.
   !---------------------------------------------------------------    
   SUBROUTINE ParticleWallContact(Particles, dt )
     
     IMPLICIT NONE
     
     TYPE(Particle_t) :: Particles
     REAL(KIND=dp) :: dt
     !---------------------------------------------------------------    
     INTEGER :: No     
     REAL(KIND=dp) :: Coord(3), Velo(3), Speed, WallVelo(3), GradVelo(3,3), Rad, Dist, &
         Mass, Normal(3), Force(3), Coeff, Spring
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: Params, BC
     TYPE(Variable_t), POINTER :: VeloVar, WallVar
     INTEGER :: i,j,k,l,n,dim, imax
     LOGICAL :: Collision,Contact,MovingWall,AnyInteraction,Found,GotVeloVar,Visited = .FALSE., &
         TrueCollision, Accumulation, AccumulationLimit, WallTrace, Stat, Hit 
     INTEGER :: Status, ElementIndex, WallNodes
     INTEGER, POINTER :: NodeIndexes(:), WallPerm(:)
     TYPE(Element_t), POINTER :: BulkElement, BoundaryElement
     TYPE(Nodes_t) :: BoundaryNodes
     REAL(KIND=dp) :: vn, dta, dtb, rn(3), v1na, v1nb, v2na, eta, AccumulationShear, ShearRate, &
         SqrtElementMetric, Weight, val, SumBasis, s
     REAL(KIND=dp), POINTER :: WallValues(:), Basis(:), dBasisdx(:,:)
     TYPE(Solver_t), POINTER :: Solver

     CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
     
     
     SAVE Visited, dim, Rad, Mass, Mesh, VeloVar, AnyInteraction, BoundaryNodes, &
         TrueCollision, WallVelo, Coeff, Velo, Coord, Spring, WallTrace, &
         WallVar, Basis, dBasisdx, GotVeloVar, Solver, Params
     
     IF( .NOT. Visited ) THEN
       Mesh => GetMesh()
       Params => GetSolverParams()
       Solver => GetSolver()

       dim = Mesh % Meshdim
       n = Mesh % MaxElementNodes
       
       ALLOCATE( Basis(n), dBasisdx(n,3) )
       
       ! Currently, one may need a different radius if the mesh leaks i.e. is 
       ! include triangles or tetrahedrans with nodes, but no faces on the surface
       !---------------------------------------------------------------    
       Rad = GetCReal( Params,'Wall Particle Radius',Found)
       IF(.NOT. Found) Rad = GetCReal( Params,'Particle Radius',Found)
       IF(.NOT. Found) THEN
         CALL Fatal('ParticleWallContact','> Particle Radius < needed!')
       END IF
       
       ! check what kind of interaction models are prescribed in the BCs
       !-----------------------------------------------------------------------
       Collision = .FALSE. 
       Contact = .FALSE.
       MovingWall = .FALSE.
       Accumulation = .FALSE.
       AccumulationLimit = .FALSE.
       WallTrace = .FALSE.

       DO k=1,CurrentModel % NumberOfBCs
         BC => CurrentModel % BCs(k) % Values
         Collision = Collision .OR. GetLogical( BC,'Wall Particle Collision',Found) 
         Contact = Contact .OR. GetLogical( BC,'Wall Particle Contact',Found) 
         MovingWall = MovingWall .OR. GetLogical( BC,'Moving Wall',Found)  
         Accumulation = Accumulation .OR. ListGetLogical( BC,'Particle Accumulation',Found)
         AccumulationLimit = AccumulationLimit .OR. ListCheckPresent( BC,'Particle Accumulation Max Shear')
         WallTrace = WallTrace .OR. ListGetLogical( BC,'Particle Trace',Found)
       END DO
       
       IF( Contact ) THEN
         Spring = GetCReal( Params,'Wall Particle Spring',Found)      
         IF(.NOT. Found) THEN
           CALL Fatal('ParticleWallContact','> Wall Particle Spring < needed!')
         END IF
       END IF
       
       IF( Collision ) THEN
         Mass = GetCReal( Params,'Particle Mass',Found)
         IF(.NOT. Found) THEN
           CALL Fatal('ParticleWallContact','> Particle Mass < needed!')
         END IF         
         Coeff = GetCReal( Params,'Wall Particle Bounciness', Found ) 
         IF(.NOT. Found) Coeff = 1.0_dp                
         TrueCollision = GetLogical( Params,'True Collision Mode',Found)      
       END IF
       

       ! Moving wall and strain based accumulation limit both require velocity. 
       ! Currently it is assumed that there can be only one velocity at a time.
       !--------------------------------------------------------------------------       
       GotVeloVar = .FALSE.
       IF( MovingWall .OR. AccumulationLimit ) THEN
         IF( MovingWall ) THEN
           IF( Collision ) THEN
             VariableName = ListGetString(Params,'Wall Velocity Variable Name',Found)
             IF( .NOT. Found ) THEN
               CALL Fatal('ParticleWallContact','Moving wall needs > Wall Velocity Variable Name <')                    
             END IF
           ELSE
             CALL Fatal('ParticleWallContact','Moving Wall assumes > Wall Particle Collision <')
           END IF
         END IF
         
         IF( AccumulationLimit ) THEN
           IF( Contact ) THEN
             VariableName = ListGetString(Params,'Velocity Variable Name',Found)
             IF( .NOT. Found ) THEN
               CALL Fatal('ParticleWallContact','Particle Accumulation needs > Velocity Variable Name <')                    
             END IF
           ELSE
             CALL Warn('ParticleWallContact','Particle Accumulation assumes > Wall Particle Collision <')
           END IF
         END IF
         
         GotVeloVar = .TRUE.
         VeloVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
         IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
           CALL Fatal('ParticleWallContact','Velocity field variable does not exist: '//TRIM(VariableName))           
         END IF
       END IF

       ! If trace is requested create the variable if it does not exist
       !---------------------------------------------------------------------
       IF( WallTrace ) THEN
         VariableName = 'Particle Trace'
         WallVar => VariableGet( Mesh % Variables,VariableName)
         IF(.NOT. ASSOCIATED(WallVar)) THEN
           ALLOCATE( WallPerm( Mesh % NumberOfNodes ) ) 
           WallPerm = 0
           CALL MakePermUsingMask(Model,Solver,Mesh,VariableName,.FALSE.,&
               WallPerm,WallNodes )

           !PRINT *,'WallNodes:',WallNodes
           IF( WallNodes > 0 ) THEN
             CALL VariableAddVector( Mesh % Variables, Mesh, Solver, VariableName, Perm = WallPerm )
           END IF
         END IF
         WallVar => VariableGet( Mesh % Variables,VariableName)
       END IF
                     
       !PRINT *,'Flags:',Contact,Collision,GotVeloVar,MovingWall, Accumulation, AccumulationLimit, WallTrace
       
       AnyInteraction = Contact .OR. Collision
       
       GradVelo = 0.0_dp
       WallVelo = 0.0_dp
       Velo = 0.0_dp
       Coord = 0.0_dp
       
       Visited = .TRUE.
     END IF
     
     IF(.NOT. AnyInteraction ) RETURN

     
     DO No = 1, Particles % NumberOfParticles

       Hit = .FALSE.
       val = 0.0_dp
       
       Status = Particles % Status(No) 
       IF( Status >= PARTICLE_LOST ) CYCLE
       IF( Status < PARTICLE_INITIATED ) CYCLE
       
       ElementIndex = Particles % ElementIndex( No ) 
       IF( ElementIndex == 0 ) CYCLE
       
       IF( Particles % InternalElements(  ElementIndex ) ) CYCLE
       
       BulkElement => Mesh % Elements( ElementIndex )
       
       IF( BulkElement % TYPE % DIMENSION == 3 ) THEN
         imax =  BulkElement % TYPE % NumberOfFaces 
       ELSE
         imax = BulkElement % TYPE % NumberOfEdges  
       END IF
       
       DO i=1, imax
         IF( BulkElement % TYPE % DIMENSION == 3 ) THEN
           j = BulkElement % FaceIndexes(i)
           BoundaryElement => Mesh % Faces( j )
         ELSE
           j = BulkElement % EdgeIndexes(i)
           BoundaryElement => Mesh % Edges(j)
         END IF
         
         Found = .FALSE.
         DO j=1,CurrentModel % NumberOfBCs
           IF(.NOT. ASSOCIATED( BoundaryElement % BoundaryInfo ) ) CYCLE
           IF ( BoundaryElement % BoundaryInfo % Constraint == &
               CurrentModel % BCs(j) % Tag ) THEN
             Found = .TRUE.
             EXIT
           END IF
         END DO
         IF( .NOT. Found ) CYCLE
         
         BC => CurrentModel % BCs(j) % Values
         Collision = GetLogical( BC,'Wall Particle Collision',Found)      
         Contact = GetLogical( BC,'Wall Particle Contact',Found)      
         IF(.NOT. (Collision .OR. Contact) ) CYCLE
         
         Coord(1:dim) = Particles % Coordinate(No,1:dim)
         CALL GetElementNodes(BoundaryNodes,BoundaryElement)
         Dist = PointFaceDistance(BoundaryElement,BoundaryNodes,Coord,Normal)
         
         ! This includes contact models using springs and possible accumulation
         ! which may be controlled by shear rate. 
         !-----------------------------------------------------------------------
         IF( Contact ) THEN
           eta = Dist - Rad
           IF( eta > 0.0 ) CYCLE
           
           Force = eta * Spring * Normal 
           Particles % Force(No,1:dim) = Particles % Force(No,1:dim) + Force(1:dim)         
           
           ! Accumulate particles that are on the boundary and then make them lost
           !-----------------------------------------------------------------------
           Accumulation = GetLogical( BC,'Particle Accumulation',Found)
           IF( Accumulation ) THEN
             stat = ParticleElementInfo( BulkElement, Coord, &
                 SqrtElementMetric, Basis, dBasisdx )
             
             AccumulationShear = GetCReal( BC,'Particle Accumulation Max Speed',Found)
             IF( Found ) THEN
               Velo(1:dim) = Particles % Velocity(No,1:dim)
               Speed = SQRT( SUM( Velo(1:dim) ** 2) )
               IF( Speed > AccumulationShear ) Accumulation = .FALSE.
             END IF
             
             AccumulationShear = GetCReal( BC,'Particle Accumulation Max Shear',Found)
             IF( Found .AND. GotVeloVar ) THEN
               CALL GetVectorFieldInMesh(VeloVar,BulkElement, Basis, WallVelo, &
                   dBasisdx, GradVelo )
               ShearRate = 0.0_dp
               DO k=1,dim
                 DO l=1,dim
                   s = 0.5 * ( GradVelo(k,l) + GradVelo(l,k) )
                   ShearRate = ShearRate + s * s
                 END DO
               END DO
               ShearRate = SQRT( ShearRate  ) 
               IF( ShearRate > AccumulationShear ) Accumulation = .FALSE.
             END IF
             
             IF( Accumulation ) THEN
               Status = PARTICLE_LOST
               Particles % Status(No) = Status
               
               Hit = .TRUE.
               WallTrace = GetLogical( BC,'Particle Trace',Found)
               val = 1.0_dp
               EXIT
             END IF
             
           END IF
           
         END IF
         
         
         ! This includes elastic and inelastic collisions that may give rise to 
         ! a contact force.
         !-----------------------------------------------------------------------
         IF( Collision ) THEN
           Velo(1:dim) = Particles % Velocity(No,1:dim)
           vn = SUM( Normal(1:dim) * Velo(1:dim) )
           
           MovingWall = GetLogical( BC,'Moving Wall',Found)
           IF( MovingWall ) THEN
             stat = ParticleElementInfo( BulkElement, Coord, &
                 SqrtElementMetric, Basis, dBasisdx )
             CALL GetVectorFieldInMesh(VeloVar,BulkElement, Basis, WallVelo )
             vn = vn - SUM( Normal(1:dim) * WallVelo(1:dim) )
           END IF
           
           IF( ABS( vn ) < TINY( vn ) ) CYCLE
           
           dta = ( Dist - Rad ) / vn
           
           IF( dta >= dt ) CYCLE        
           
           ! A historical collision but the distance is growing
           IF( dta < 0.0 .AND. Dist - Rad > 0.0 ) CYCLE      
           
           IF( dta < 0.0 ) THEN
             PRINT *,'Coord:',Coord(1:dim)
             PRINT *,'Velo:',Velo(1:dim)
             PRINT *,'Dist',Dist,Rad
             PRINT *,'vn',vn,dta,dt
           END IF
           
           ! these are defined as is so that we could reuse the binary particle collision stuff
           ! at the limit Mass2 -> infinity
           dtb = dt - dta
           rn = -Normal
           v1na = SUM( Velo(1:dim) * rn(1:dim) )
           v2na = SUM( WallVelo(1:dim) * rn(1:dim) ) 
           
           v1nb = Coeff * (v2na - v1na) + v2na 
           
           ! Set either force or velocity directly 
           ! only the normal component of velocity/force is affected by collisions    
           IF( TrueCollision ) THEN
             ! compute the path until the collision
             Coord = Coord + dta * Velo
             Velo = Velo + (v1nb-v1na) * rn
             ! compute the path after the collision
             Coord = Coord + dtb * Velo
           ELSE
             Coord = Coord + (v1na-v1nb) * rn * dta
             Force = Mass * (v1nb-v1na) * rn / dt
           END IF
           
           IF( dta < 0.0 ) THEN
             PRINT *,'dtb',dtb
             PRINT *,'rn',rn
             PRINT *,'na',v1na,v2na,v1nb
             PRINT *,'Coord2',Coord
             PRINT *,'Force2',Force
           END IF
           
           IF( TrueCollision ) THEN
             Particles % Coordinate(No,1:dim) = Coord(1:dim)
             Particles % Velocity(No,1:dim) = Velo(1:dim)
           ELSE
             Particles % Coordinate(No,1:dim) = Coord(1:dim)
             Particles % Force(No,1:dim) = Particles % Force(No,1:dim) + Force(1:dim)         
           END IF
           
           Hit = .TRUE.
           WallTrace = GetLogical( Params,'Particle Trace',Found)
           IF( WallTrace) THEN
             IF( .NOT. MovingWall ) THEN
               ! we need the local basis functions but they haven't been computed
               stat = ParticleElementInfo( BulkElement, Coord, &
                   SqrtElementMetric, Basis, dBasisdx )
             END IF
             val = Mass * (v1nb - v1na ) / dt
           END IF
           
           ! Only one collision for each particle & element
           EXIT
         END IF
       END DO


       IF(.NOT. Hit) CYCLE

       ! If requested populate a result vector that shows the 
       ! total hits or forces on the boundary.
       !------------------------------------------------------
       IF( WallTrace ) THEN
         n = BulkElement % TYPE % NumberOfNodes
         NodeIndexes => BulkElement % NodeIndexes
         SumBasis = 0.0_dp
         DO i = 1, n
           j = WallVar % Perm( NodeIndexes(i) )
           IF( j > 0 ) SumBasis = SumBasis + Basis(i)
         END DO
         
         ! only populate the vector if there really where some 
         ! hits resulting to nonzero sum of the basis vectors.
         ! The sum at boundary is normalized to one.
         !----------------------------------------------------
         IF( SumBasis > TINY( SumBasis ) ) THEN
           val = val / SumBasis
           DO i = 1,n
             j = WallVar % Perm( NodeIndexes(i) )
             IF( j > 0 ) THEN
               weight = Basis(i)
               WallVar % Values( j ) = WallVar % Values( j ) + weight * val
             END IF
           END DO
         END IF
       END IF
     END DO

   END SUBROUTINE ParticleWallContact


!------------------------------------------------------------------------------
END SUBROUTINE ParticleDynamics
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialize the particle dynamics solver with some defaults.
!------------------------------------------------------------------------------
SUBROUTINE ParticleDynamics_Init( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE Lists

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------
  
  TYPE(ValueList_t), POINTER :: Params
  LOGICAL :: Found

  Params => Solver % Values

END SUBROUTINE ParticleDynamics_Init

!> \}

