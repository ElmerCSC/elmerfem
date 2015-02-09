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

!-------------------------------------------------------------------------------
!> Subroutine for advecting fields in time using particles to follow them  
!> backwards in time, and taking the field value from the given point. This should overcome
!> all problems with diffusion. 
!> This is a dynamically loaded solver with a standard interface. 
!> \ingroup Solvers
!-------------------------------------------------------------------------------
SUBROUTINE ParticleAdvector( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils

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
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Solver_t), POINTER :: PSolver
  TYPE(Variable_t), POINTER :: Var, PTimeVar
  LOGICAL :: GotIt, Debug, Hit, InitLocation, InitTimestep, Found, ParticleInfo
  INTEGER :: i,j,k,n,dim,No,nodims,&
      ElementIndex, VisitedTimes = 0, nstep, &
      Status,TimeOrder, PartitionChanges, TimeStepsTaken=0,&
      ParticleStepsTaken=0, TotParticleStepsTaken, TotNoParticles, &
      istep,iorder
  REAL(KIND=dp) :: maxdt, dertime = 0.0, tottime = 0.0
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, IntegMethod
  TYPE(Particle_t), POINTER  :: Particles
  
  SAVE Nstep, VisitedTimes, TimeOrder, &
      tottime, TimeStepsTaken, ParticleStepsTaken, ParticleInfo

!------------------------------------------------------------------------------

  CALL Info('ParticleAdvector','-----------------------------------------', Level=4 )
  CALL Info('ParticleAdvector','Advecting fields using particle tracking',Level=4) 
  
  Particles => GlobalParticles
  VisitedTimes = VisitedTimes + 1

  Params => Solver % Values
  Mesh => Solver % Mesh
  PSolver => Solver
  DIM = CoordinateSystemDimension()

  maxdt = 0.0_dp  
  istep = 1
  iorder = 1

  ! Do some initalialization: allocate space, check fields  
  !------------------------------------------------------------------------
  IF( VisitedTimes == 1 ) THEN
    TimeOrder = GetInteger( Params,'Time Order',GotIt)
    CALL SetParticlePreliminaries( Particles, dim, TimeOrder )
    Nstep = GetInteger( Params,'Max Timestep Intervals',Found)
    IF(.NOT. Found) Nstep = 1
    ParticleInfo = GetLogical( Params,'Particle Info',Found)
  END IF

  ! Initialize particles always since we just advance -dt each time
  !-------------------------------------------------------------------------
  IF( VisitedTimes == 1 .OR. GetLogical( Params,'Reinitialize Particles',GotIt ) ) THEN
    CALL InitializeParticles( Particles ) 
    CALL ReleaseWaitingParticles(Particles) 
    Particles % Status = PARTICLE_LOCATED
  ELSE
    ! in case the velocity field is changed update also the particle velocities
    CALL SetParticleVelocities()
  END IF

  IF( VisitedTimes == 1 ) THEN
    IF( GetLogical( Params,'Particle Time',Found) ) THEN
      CALL ParticleVariableCreate( Particles,'particle time')
    END IF 
    CALL ParticleVariableCreate( Particles,'particle distance')
  ELSE	 
    PtimeVar => ParticleVariableGet( Particles, 'particle time' )
    IF( ASSOCIATED( PTimeVar ) ) THEN
      PTimeVar % Values = 0.0_dp
    END IF
  END IF

  ! Freeze particles that are known not to move (e.g. no-slip wall)
  !----------------------------------------------------------------
  CALL SetFixedParticles( )

  ! We are advecting backwards in time!
  !----------------------------------------------------
  InitTimestep = .TRUE.
  IF( Particles % RK2 ) THEN
    iorder = 2
  ELSE
    iorder = 1
  END IF

!  CALL ParticleStatusCount( Particles )


  DO i=1,nstep

    ! Get the timestep size, initialize at 1st round
    !--------------------------------------------------------------
    maxdt = GetParticleTimeStep( Particles, InitTimestep )

    ! If size of timestep goes to zero then no more steps are needed
    !---------------------------------------------------------------
    IF( ABS( maxdt ) < TINY( maxdt ) ) THEN
      WRITE (Message,'(A,I0)') 'Number of steps used: ',i-1
      CALL Info('ParticleAdvector',Message,Level=6)
      EXIT	
    END IF

    dertime = dertime + maxdt
    tottime = tottime + maxdt

    TimeStepsTaken = TimeStepsTaken + 1
    ParticleStepsTaken = ParticleStepsTaken + Particles % NumberOfParticles

    ! If there are periodic BCs apply them just before locating the particles
    !------------------------------------------------------------------------
    ! CALL ParticleBoxPeriodic( Particles )


    ! Loop over possible runge-kutta steps
    !------------------------------------------------------------------------

    DO istep = 1, iorder

      ! Set velocity field at points. This is not done the first time. 
      ! This has to be here since the velocity field could have changed
      ! between calls. 
      !------------------------------------------------------------------
      !IF( .NOT. InitTimestep ) CALL SetParticleVelocities()

      ! Advance the coordinate, r = r0 + v * dt
      ! either with 2nd order Runge-Kutta, or with 1st order explicit scheme.
      !------------------------------------------------------------------
      CALL ParticleAdvanceTimestep( Particles, istep )

      !    CALL ParticleStatusCount( Particles )

      ! Find the elements (and only the elements) in which the particles are in. 
      !------------------------------------------------------------------------    
      CALL LocateParticles( Particles ) 

      CALL SetParticleVelocities()

      ! Integrate over the particle path (\int f(r) ds or \int f(r) dt )
      !------------------------------------------------------------------
      CALL ParticlePathIntegral( Particles, istep )

      InitTimestep = .FALSE.
    END DO 

    WRITE (Message,'(A,I0,A,I0,A)') 'Timestep ',i,' with ',&
	Particles % NumberOfMovingParticles,' moving particles'
    CALL Info('ParticleAdvector',Message,Level=6)

    !CALL ParticleInformation(Particles, ParticleStepsTaken, &
    !	TimeStepsTaken, tottime )

  END DO

  ! Set the advected field giving the final locations of the particles backward in time
  !------------------------------------------------------------------------------------
  CALL SetAdvectedField()
  
  ! In the end show some statistical info
  !---------------------------------------------------------------   
  IF( ParticleInfo ) THEN
    CALL ParticleInformation(Particles, ParticleStepsTaken, &
	TimeStepsTaken, tottime )
  END IF
    
  CALL Info('ParticleAdvector','All done',Level=4)
  CALL Info('ParticleAdvector', '-----------------------------------------', Level=4 )
  

CONTAINS

  !> Eliminate particles that sit an fixed boundaries.
  !-------------------------------------------------------------------------
  SUBROUTINE SetFixedParticles( )

    INTEGER :: i,j,k
    TYPE( ValueList_t), POINTER :: BC, BodyForce
    TYPE( Element_t ), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: PTime(:), PCond(:)
    REAL(KIND=dp) :: PTimeConst
    LOGICAL :: Found, GotTime, SomeBodyForce, SomeBC
    TYPE(ValueList_t), POINTER :: Params


    Params => GetSolverParams()

    SomeBC = ListCheckPresentAnyBC( Model,'Particle Fixed Condition') 
    SomeBodyForce = ListCheckPresentAnyBodyForce( Model,'Particle Fixed Condition') 

    IF( .NOT. (SomeBC .OR. SomeBodyForce ) ) RETURN		  

    i = Model % MaxElementNodes 
    ALLOCATE( PTime(i), PCond(i) )

    PtimeVar => ParticleVariableGet( Particles, 'particle time' )
    GotTime = ASSOCIATED( PTimeVar ) 
    PTimeConst = ListGetCReal( Params,'Fixed Particle Time',Found)

    IF( SomeBC ) THEN
      DO j=Mesh % NumberOfBulkElements+1,&
  	Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(j)
        Model % CurrentElement => Element
        BC => GetBC( Element )
	n = GetElementNOFNodes()

        PCond(1:n) = GetReal( BC,'Particle Fixed Condition',Found)
        IF(.NOT. Found ) CYCLE
	
        PTime(1:n) = GetReal( BC,'Particle Time',Found)
    
        DO i=1,n
          IF( PCond(i) < 0.0_dp ) CYCLE

	  k = Element % NodeIndexes(i)
          Particles % Status( k ) = PARTICLE_FIXEDCOORD
	  IF(.NOT. GotTime ) CYCLE

	  IF( Found ) THEN
            PTimeVar % Values( k ) = PTime(i)
          ELSE
            PTimeVar % Values( k ) = PTimeConst
          END IF 
        END DO
      END DO

      i = COUNT( Particles % Status == PARTICLE_FIXEDCOORD )
      WRITE( Message,'(A,I0)') 'Number of fixed boundary particles: ',i 
      CALL Info('SetFixedParticles',Message)
    END IF 


    IF( SomeBodyForce ) THEN
      DO j=1, Mesh % NumberOfBulkElements 
        Element => Mesh % Elements(j)
        Model % CurrentElement => Element
        BodyForce => GetBodyForce( Element )
	n = GetElementNOFNodes()

        PCond(1:n) = GetReal( BodyForce,'Particle Fixed Condition',Found)
        IF(.NOT. Found ) CYCLE
	
        PTime(1:n) = GetReal( BodyForce,'Particle Time',Found)
    
        DO i=1,n
          IF( PCond(i) < 0.0_dp ) CYCLE

	  k = Element % NodeIndexes(i)
          Particles % Status( k ) = PARTICLE_FIXEDCOORD
	  IF(.NOT. GotTime ) CYCLE

	  IF( Found ) THEN
            PTimeVar % Values( k ) = PTime(i)
          ELSE
            PTimeVar % Values( k ) = PTimeConst
          END IF 
        END DO
      END DO

      i = COUNT( Particles % Status == PARTICLE_FIXEDCOORD )
      WRITE( Message,'(A,I0)') 'Number of fixed bulk particles: ',i 
      CALL Info('SetFixedParticles',Message)
    END IF 

   DEALLOCATE( PTime, PCond )

  END SUBROUTINE SetFixedParticles 

  
  !------------------------------------------------------------------------
  !> Compute field values at the given points in the FE mesh. 
  !-------------------------------------------------------------------------
  SUBROUTINE SetParticleVelocities()

    TYPE(Element_t), POINTER :: BulkElement
    INTEGER :: No, Status
    REAL(KIND=dp) :: Coord(3),Velo(3),GradVelo(3,3)    
    TYPE(Element_t), POINTER :: BulkElement2
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: Params
    REAL(KIND=dp) :: VeloAtPoint(3), GradVeloAtPoint(3,3),dtime
    LOGICAL :: Stat, UseGradVelo, Visited = .FALSE.
    INTEGER :: i,j,k,l,n,dim,TimeOrder, NewLost(3), OldLost, FixedLost
    INTEGER, POINTER :: NodeIndexes(:), FieldPerm(:),FieldIndexes(:)
    REAL(KIND=dp) :: SqrtElementMetric, Weight, Speed, SpeedMin
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:), Coordinate(:,:), Velocity(:,:)
    LOGICAL :: GotIt
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    TYPE(Variable_t), POINTER :: VeloVar
    TYPE(Variable_t), POINTER :: DtVar	
    
    
    SAVE :: Visited, Mesh, dim, Basis, dBasisdx, Params, VeloVar, UseGradvelo, DtVar, &
	SpeedMin, NewLost 

    IF( .NOT. Visited ) THEN
      Mesh => GetMesh()
      dim = Mesh % MeshDim
      n = Mesh % MaxElementNodes
      ALLOCATE( Basis(n), dBasisdx(n, 3) )
      
      Params => GetSolverParams()
      
      VariableName = ListGetString(Params,'Velocity Variable Name',Found)
      IF(.NOT. Found) VariableName = 'flow solution'
      VeloVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
      IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
        CALL Fatal('ParticleFieldInteraction','Velocity field variable does not exist: '//TRIM(VariableName))           
      END IF
      UseGradVelo = GetLogical( Params,'Velocity Gradient Correction',Found)

      IF(UseGradVelo .AND. Particles % RK2 ) THEN
        CALL Warn('ParticlePathIntegral','Quadratic source correction incompatibe with Runge-Kutta')
        UseGradVelo = .FALSE.
      END IF

      IF( .NOT. Particles % DtConstant ) THEN
        DtVar => ParticleVariableGet( Particles,'particle dt')
        IF( .NOT. ASSOCIATED( DtVar ) ) THEN
          CALL Fatal('SetParticleVelocities','Required field > particle dt < not present!')
        END IF
      END IF      

      SpeedMin = ListGetConstReal( Params,'Particle Min Speed',Found)
      IF(.NOT. Found) SpeedMin = EPSILON( SpeedMin )

      NewLost =  0

      Visited = .TRUE.
    END IF 

    Coordinate => Particles % Coordinate
    Velocity => Particles % Velocity
    Coord = 0.0_dp
    Velo = 0.0_dp

    OldLost = 0
    FixedLost = 0

    IF( Particles % DtConstant ) THEN
      dtime = Particles % DtSign * Particles % dtime
    END IF

    DO No = 1, Particles % NumberOfParticles
      Status = GetParticleStatus( Particles, No )
      IF( Status >= PARTICLE_LOST .OR. &
          Status <= PARTICLE_INITIATED .OR. &
          Status == PARTICLE_FIXEDCOORD .OR. &
          Status == PARTICLE_WALLBOUNDARY ) THEN
        OldLost = OldLost + 1
	CYCLE
      END IF
      
      ElementIndex = GetParticleElement( Particles, No )
      IF( ElementIndex == 0 ) THEN
!PRINT *,'elemindex  = 0'
        Particles % Status(No) = PARTICLE_LOST
        NewLost(1) = NewLost(1) + 1
        CYCLE       
      END IF

      BulkElement => Mesh % Elements( ElementIndex )
      
      Coord(1:dim) = Coordinate( No, 1:dim )

      !-------------------------------------------------------------------------
      ! Get velocity from mesh
      !-------------------------------------------------------------------------
      IF( UseGradVelo ) THEN
        stat = ParticleElementInfo( BulkElement, Coord, &
            SqrtElementMetric, Basis, dBasisdx )
      ELSE
        stat = ParticleElementInfo( BulkElement, Coord, &
            SqrtElementMetric, Basis )
      END IF

      IF(.NOT. Stat ) THEN
!	print *,'Particle not in element!',No,Coord
        Particles % Status(No) = PARTICLE_LOST
        NewLost(2) = NewLost(2) + 1
        CYCLE
      END IF

      IF( UseGradVelo ) THEN
        CALL GetVectorFieldInMesh(VeloVar,BulkElement, Basis, VeloAtPoint, &
            dBasisdx, GradVeloAtPoint )
	IF( .NOT. Particles % DtConstant ) THEN
          dtime = Particles % DtSign * DtVar % Values(No)
        END IF
        DO i=1,dim
          Velo(i) = VeloAtPoint(i) + &
              0.5_dp * SUM( GradVeloAtPoint(i,1:dim) * VeloAtPoint(1:dim) ) * dtime        
        END DO
      ELSE
        CALL GetVectorFieldInMesh(VeloVar, BulkElement, Basis, VeloAtPoint )
	Velo(1:dim) = VeloAtPoint(1:dim)
      END IF

      Speed = SQRT( SUM( Velo(1:dim) ** 2 ) )
      IF( Speed < SpeedMin ) THEN
 	Particles % Status(No) = PARTICLE_FIXEDCOORD
        Velocity( No, 1:dim ) = 0.0_dp
        FixedLost = FixedLost + 1
      ELSE
        Velocity( No, 1:dim ) = Velo(1:dim)
      END IF
    END DO

    IF( .FALSE. ) THEN
      IF( NewLost(1) > 0 ) CALL Warn('SetParticleVelocities','Some particles could not be located')
      PRINT *,'Total number of particles:',Particles % NumberOfParticles
      PRINT *,'Passive particles:',OldLost
      PRINT *,'New lost particles:',NewLost
      PRINT *,'New fixed velo particles:',FixedLost
    END IF

!    print *,'NewLost:',NewLost

    
  END SUBROUTINE SetParticleVelocities
   



  !------------------------------------------------------------------------
  !> Compute field values at the given points in the FE mesh. 
  !-------------------------------------------------------------------------
  SUBROUTINE SetAdvectedField()

    TYPE(Element_t), POINTER :: BulkElement
    INTEGER :: No, Status, NoParticles
    REAL(KIND=dp) :: dtime, Coord(3),Velo(3),val,vals(10)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: Params
    LOGICAL :: Stat, Visited = .FALSE.
    INTEGER :: i,j,k,l,n,nsize,dim,wallcount,NoVar,NoNorm,dofs,maxdim
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: SqrtElementMetric, Norm, PrevNorm = 0.0_dp, Change
    REAL(KIND=dp), POINTER :: Basis(:)
    LOGICAL :: GotIt, Difference,Cumulative,Derivative,GotVar,GotRes,GotOper,Debug,&
        UsePerm,InternalVariable,Initiated, Parallel
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, ResultName, OperName, Name
    TYPE(Variable_t), POINTER :: TargetVar, ResultVar, Var
    TYPE(Variable_t), POINTER :: ParticleVar
    REAL(KIND=dp), POINTER :: TmpValues(:), NodeValues(:), NewValues(:)
    INTEGER, POINTER :: TmpPerm(:), UnitPerm(:)
    REAL(KIND=dp) :: x,y,z
    
    
    SAVE :: Visited, PrevNorm, UnitPerm
    
    Mesh => GetMesh()
    dim = Mesh % MeshDim
    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )
    Coord = 0.0_dp
    Velo = 0.0_dp

    Params => GetSolverParams()
    NoNorm = GetInteger( Params,'Norm Variable Index',GotIt)
    NoParticles = Particles % NumberOfParticles
    maxdim = 0

    Parallel = ( ParEnv % PEs > 1 ) 

    Initiated = .FALSE.
100 NoVar = 0 
    DO WHILE(.TRUE.)
      NoVar = NoVar + 1

      WRITE (Name,'(A,I0)') 'Variable ',NoVar
      VariableName = GetString( Params,Name,GotVar)
      IF(.NOT. GotVar ) EXIT

      ! Get the target variables
      ! Variables starting with 'particle' as associated with particles
      !----------------------------------------------------------------
      IF( VariableName == 'particle coordinate' .OR. &
          VariableName == 'particle velocity' .OR. &
          VariableName == 'particle force') THEN
        dofs = dim 
        InternalVariable = .TRUE.
        maxdim = MAX( dim, maxdim )
      ELSE IF( VariableName(1:8) == 'particle' ) THEN
        dofs = 1
        InternalVariable = .TRUE.
        maxdim = MAX( 1, maxdim )
      ELSE
        TargetVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
        IF( .NOT. ASSOCIATED(TargetVar)) THEN
          CALL Fatal('ParticleAdvector','Could not obtain target variable:'&
              //TRIM(VariableName))
          CYCLE
        END IF
        dofs = TargetVar % dofs
        IF( dofs > 1 ) THEN
          CALL Fatal('ParticleAdvector','Advection implemented so far only for scalars')
        END IF
        InternalVariable = .FALSE.
        maxdim = MAX( dofs, maxdim ) 
      END IF

      IF(.NOT. Initiated ) CYCLE


      ! For internal variables the target name is the field name
      !---------------------------------------------------------      
      Difference = .FALSE.
      Derivative = .FALSE.
      Cumulative = .FALSE.

      WRITE (Name,'(A,I0)') 'Operator ',NoVar
      OperName = GetString( Params,Name,GotOper)
      IF( GotOper ) THEN
        IF( OperName == 'difference' ) THEN
          Difference = .TRUE.
        ELSE IF( OperName == 'derivative') THEN
          Derivative = .TRUE.
        ELSE IF( OperName == 'cumulative' .OR. OperName == 'age' ) THEN
          Cumulative = .TRUE.
        ELSE
          CALL Fatal('SetAdvectedField','Unknown operator: '//TRIM(OperName))
        END IF
      ELSE
        OperName = 'adv'
      END IF
      
      WRITE (Name,'(A,I0)') 'Result Variable ',NoVar
      ResultName = GetString( Params,Name,GotRes)
      IF( .NOT. GotRes ) THEN
        ResultName = TRIM(OperName)//' '//TRIM(VariableName)
      END IF

      ! Create variables if they do not exist
      !---------------------------------------------------------      
      ResultVar => VariableGet( Mesh % Variables, TRIM(ResultName) )
      IF( .NOT. ASSOCIATED(ResultVar)) THEN
        IF( InternalVariable ) THEN
          UsePerm = .FALSE.
        ELSE
          ! Inherit the Perm from the target variable
          UsePerm = ASSOCIATED( TargetVar % Perm ) 
        END IF

        ! This trick is done to allow postprocessing routines to work better
        IF(.NOT. UsePerm ) THEN
          IF(.NOT. ASSOCIATED(UnitPerm)) THEN
            ALLOCATE( UnitPerm(Mesh % NumberOfNodes ) )
            DO i=1,Mesh % NumberOfNodes
              UnitPerm(i) = i
            END DO
          END IF
        END IF

        IF(UsePerm) THEN
          CALL VariableAddVector( Mesh % Variables,Mesh,PSolver,ResultName,dofs,&
              Perm = TargetVar % Perm )
        ELSE
          CALL VariableAddVector( Mesh % Variables,Mesh,PSolver,ResultName,dofs, &
              Perm = UnitPerm )
        END IF
        IF( dofs == 1 ) THEN
          CALL Info('ParticleAdvector','Created a scalar variable: '//TRIM(ResultName) )
        ELSE
          CALL Info('ParticleAdvector','Created a vector variable: '//TRIM(ResultName) )          
        END IF
        ResultVar => VariableGet( Mesh % Variables, TRIM(ResultName))
        IF(.NOT. ASSOCIATED(ResultVar)) CALL Fatal('ParticleAdvector','Problems in VariableAdd')
      END IF
      

      ! Finally, set the values
      !---------------------------------------------------------      
      IF( InternalVariable ) THEN
        IF( VariableName == 'particle coordinate') THEN
          IF( ResultVar % Dofs /= dim ) THEN
            CALL Fatal('ParticleAdvector','Variable should have dim dofs: '//TRIM(VariableName))
          END IF
          DO i=1,NoParticles
            DO j=1,dim
              NewValues(dim*(i-1)+j) = Particles % Coordinate(i,j)
            END DO
          END DO
        ELSE IF( VariableName == 'particle coordinate_abs') THEN
          DO i=1,NoParticles
            val = 0.0_dp
            DO j=1,dim
              val = val + Particles % Coordinate(i,j)**2
            END DO
            NewValues(i) = SQRT( val )
          END DO
        ELSE IF( VariableName == 'particle velocity') THEN
          IF( ResultVar % Dofs /= dim ) THEN
            CALL Fatal('ParticleAdvector','Variable should have dim dofs: '//TRIM(VariableName))
          END IF
          DO i=1,NoParticles
            DO j=1,dim
              NewValues(dim*(i-1)+j) = Particles % Velocity(i,j)
            END DO
          END DO
          
        ELSE IF( VariableName == 'particle velocity_abs') THEN
          DO i=1,NoParticles
            val = 0.0_dp
            DO j=1,dim
              val = val + Particles % Velocity(i,j)**2
            END DO
            NewValues(i) = SQRT( val )
          END DO

        ELSE IF( VariableName == 'particle force') THEN
          IF( ResultVar % Dofs /= dim ) THEN
            CALL Fatal('ParticleAdvector','Variable should have dim dofs: '//TRIM(VariableName))
          END IF
          DO i=1,NoParticles
            DO j=1,dim
              NewValues(dim*(i-1)+j) = Particles % Force(i,j)
            END DO
          END DO

        ELSE IF( VariableName == 'particle status') THEN
          DO i=1,NoParticles
            NewValues(i) = 1.0_dp * Particles % Status(i)
          END DO
          
        ELSE IF( VariableName(1:8) == 'particle' ) THEN
          ParticleVar => ParticleVariableGet( Particles, VariableName )
          IF( ASSOCIATED( ParticleVar ) ) THEN
            NewValues = ParticleVar % Values
          ELSE
            CALL Warn('ParticleAdvector','Field does not exist: '//TRIM(VariableName))
          END IF
        END IF
        
      ELSE 

        DO i = 1, NoParticles
          Status = GetParticleStatus( Particles, i )
          
          IF( Status >= PARTICLE_LOST ) CYCLE
          IF( Status <= PARTICLE_INITIATED ) CYCLE
          
          ElementIndex = GetParticleElement( Particles, i )

          IF( ElementIndex == 0 ) CYCLE

          BulkElement => Mesh % Elements( ElementIndex )      
          
          Coord(1:dim) = Particles % Coordinate(i, 1:dim) 
          Velo(1:dim) = Particles % Velocity(i, 1:dim) 
          
          stat = ParticleElementInfo( BulkElement, Coord, SqrtElementMetric, Basis )
          IF(.NOT. stat) CYCLE
          
          IF( dofs == 1 ) THEN
            CALL GetScalarFieldInMesh(TargetVar, BulkElement, Basis, val ) 
            NewValues( i ) =  val 
          ELSE
            CALL GetVectorFieldInMesh(TargetVar, BulkElement, Basis, vals ) 
            DO j=1,dofs             
              NewValues( dofs*(i-1)+j ) = vals(j)     
            END DO
          END IF
        END DO
      END IF

      ! In a serial case the nodes and particles are directly associated. 
      ! In a parallel case we need to transfer the values from particles in 
      ! different partitions to nodes. 
      !---------------------------------------------------------------------
      IF( Parallel ) THEN
        CALL ParticleAdvectParallel( Particles, NewValues, NodeValues, dofs )
      END IF

      ! Finally move the nodal values to the target variable 
      !---------------------------------------------------------------------
      DO j=1,Mesh % NumberOfNodes 
        k = j
        IF( ASSOCIATED( ResultVar % Perm ) ) k = ResultVar % Perm( k )
        IF( k == 0 ) CYCLE
        IF( Difference .OR. Derivative ) THEN
          ResultVar % Values( k ) = NodeValues( j ) - TargetVar % Values( k ) 
        ELSE IF( Cumulative ) THEN
          ResultVar % Values( k ) = NodeValues( j ) + ResultVar % Values( k )
        ELSE
          ResultVar % Values( k ) = NodeValues( j ) 
        END IF
      END DO
      IF( Derivative ) ResultVar % Values = ResultVar % Values / dertime 

      ! To allow computation of change in the standard manner the Variable
      ! is set to point to the one of interest. This is mainly used in the 
      ! tests, or could be used in for convergence monitoring also. 
      !---------------------------------------------------------------
      IF( NoVar == NoNorm ) THEN
        n = SIZE( ResultVar % Values ) 
!        Norm = ComputeNorm( Solver, n, ResultVar % Values ) 
        Norm = SQRT( SUM( ResultVar % Values ** 2) / n )
        Change = 2.0 * ABS( Norm-PrevNorm ) / ( Norm + PrevNorm )
        PrevNorm = Norm

        Solver % Variable % Norm = Norm
        Solver % Variable % NonlinChange = Change

        ! Here the name is ComputeChange in order to get the change also to ElmerGUI
        ! albeit in a very dirt style. One could also edit ElmerGUI....
        WRITE( Message, '(a,g15.8,g15.8,a)') &
            'SS (ITER=1) (NRM,RELC): (',Norm, Change,' ) :: '//TRIM( ResultName )
        CALL Info( 'ComputeChange', Message, Level=3 )
      END IF
    END DO


    ! Allocate the local new temporal values
    IF(.NOT. Initiated ) THEN
      NoVar = NoVar - 1
      IF( NoVar < 1 ) THEN
        CALL Fatal('ParticleAdvector','No target and result variables exist!')
      END IF
      ALLOCATE( NewValues( maxdim * Particles % NumberOfParticles ) ) 
      NewValues = 0.0_dp
      IF( Parallel ) THEN
        ALLOCATE( NodeValues( maxdim * Mesh % NumberOfNodes ) )
        NodeValues = 0.0_dp
      ELSE
        NodeValues => NewValues
      END IF
      Initiated = .TRUE.
      GOTO 100
    END IF

    DEALLOCATE( NewValues ) 
    IF( Parallel ) DEALLOCATE( NodeValues ) 
    Visited = .TRUE.

  END SUBROUTINE SetAdvectedField

!------------------------------------------------------------------------------
END SUBROUTINE ParticleAdvector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Initialization for the primary solver: ParticleAdvector
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE ParticleAdvector_Init( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE Interpolation
  USE MeshUtils
  USE ElementUtils
  USE ParticleUtils

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
  INTEGER :: NormInd

  Params => Solver % Values

  ! These are default setting that make the operation of the advection solver 
  ! possible. There should always be one passive particle for each active node.
  !---------------------------------------------------------------------------
  CALL ListAddString( Params,'Coordinate Initialization Method','nodal ordered')
  CALL ListAddString( Params,'Velocity Initialization Method','nodal velocity')
  CALL ListAddInteger( Params,'Time Order',0 )
  CALL ListAddConstReal( Params,'Particle Node Fraction',1.0_dp)
  IF(.NOT. ListCheckPresent( Params,'Particle Accurate At Face') ) &
      CALL ListAddLogical( Params,'Particle Accurate At Face',.TRUE.)  
  CALL ListAddLogical( Params,'Particle Dt Negative',.TRUE.)
  CALL ListAddLogical( Params,'Particle Fix Frozen',.TRUE.)

  ! If we want to show a pseudonorm add a variable for which the norm
  ! is associated with.
  NormInd = ListGetInteger( Params,'Show Norm Index',Found)
  IF( NormInd > 0 ) THEN
    IF( .NOT. ListCheckPresent( Params,'Variable') ) THEN
      CALL ListAddString( Solver % Values,'Variable','-nooutput -global particleadvector_var')
    END IF
  END IF


END SUBROUTINE ParticleAdvector_Init
