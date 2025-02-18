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
! *  Authors: Peter Raback & Juha Ruokolainen
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
  LOGICAL :: GotIt, Hit, InitLocation, InitTimestep, Found, &
      ParticleInfo, InitAllVelo, ReverseTime, Reinitialize
  INTEGER :: i,j,k,n,dim,No,nodims,&
      ElementIndex, VisitedTimes = 0, nstep, &
      Status,TimeOrder, PartitionChanges, TimeStepsTaken=0,&
      ParticleStepsTaken=0, TotParticleStepsTaken, TotNoParticles, &
      istep,iorder,NoMoving
  REAL(KIND=dp) :: maxdt, dertime = 0.0, tottime = 0.0
  REAL(KIND=dp), ALLOCATABLE :: OrigCoordinate(:,:)
  CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, IntegMethod
  TYPE(Particle_t), POINTER  :: Particles
  CHARACTER(*), PARAMETER :: Caller = 'ParticleAdvector'

  
  SAVE Nstep, VisitedTimes, TimeOrder, &
      tottime, TimeStepsTaken, ParticleStepsTaken, ParticleInfo

!------------------------------------------------------------------------------

  CALL Info(Caller,'-----------------------------------------', Level=4 )
  CALL Info(Caller,'Advecting fields using particle tracking',Level=4) 
  
  Particles => GlobalParticles
  VisitedTimes = VisitedTimes + 1

  Params => GetSolverParams()
  Mesh => Solver % Mesh
  PSolver => Solver
  DIM = CoordinateSystemDimension()

  ReverseTime = ListGetLogical( Params,'Particle time reverse',GotIt)
  Reinitialize = GetLogical( Params,'Reinitialize Particles',GotIt ) 
  IF( ReverseTime .AND. .NOT. Reinitialize) THEN
    IF(GotIt) THEN
      CALL Fatal(Caller,'Reverse time and no reinitialization are in conflict!')
    ELSE
      CALL Info(Caller,'Enforcing reinitialization because of reversed time!')
      Reinitialize = .TRUE.
    END IF
  END IF

  IF( ListGetLogical( Params,'Fix Tangent Velocity',GotIt ) ) THEN
    CALL ApplyTangentFix(.FALSE.)
  END IF
  
  
  maxdt = 0.0_dp  
  istep = 1
  iorder = 1

  InitAllVelo = .TRUE.
  
  ! Do some initialization: allocate space, check fields
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
  IF( VisitedTimes == 1 .OR. Reinitialize ) THEN
    CALL InitializeParticles( Particles, IsAdvector = .TRUE.) 
    CALL ReleaseWaitingParticles(Particles) 
    Particles % Status = PARTICLE_LOCATED
    
    IF( ReverseTime ) THEN
      n = Particles % NumberOfParticles
      ALLOCATE(OrigCoordinate(n,dim))      
      DO i=1,n
        DO j=1,dim
          OrigCoordinate(i,j) = Particles % Coordinate(i,j) 
        END DO
      END DO
    END IF
  ELSE
    ! in case the velocity field is changed update also the particle velocities
    CALL SetParticleVelocities(InitAllVelo)
    InitAllVelo = .FALSE.
  END IF

  
  IF( VisitedTimes == 1 ) THEN
    CALL ParticleVariableCreate( Particles,'particle time')
    CALL ParticleVariableCreate( Particles,'particle distance')
  END IF
  PtimeVar => ParticleVariableGet( Particles, 'particle time' )
  IF( ASSOCIATED( PTimeVar ) ) THEN
    PTimeVar % Values = 0.0_dp
  END IF
  Var => ParticleVariableGet( Particles, 'particle distance' )
  IF( ASSOCIATED( Var ) ) THEN
    Var % Values = 0.0_dp
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

  ! This is the default for particle advector!
  Particles % DtSign = -1

1 CONTINUE
  
  DO i=1,nstep
    ! Get the timestep size, initialize at 1st round
    !--------------------------------------------------------------
    maxdt = GetParticleTimeStep( Particles, InitTimestep )

    ! If size of timestep goes to zero then no more steps are needed
    !---------------------------------------------------------------
    IF( ABS( maxdt ) < TINY( maxdt ) ) THEN
      WRITE (Message,'(A,I0)') 'Number of steps used: ',i-1
      CALL Info(Caller,Message,Level=6)
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

      IF( InfoActive( 20 ) ) THEN
        CALL ParticleStatusCount( Particles )
      END IF
        
      ! Find the elements (and only the elements) in which the particles are in. 
      !------------------------------------------------------------------------    
      CALL LocateParticles( Particles ) 

      CALL SetParticleVelocities(InitAllVelo)
      InitAllVelo = .FALSE.

      ! Integrate over the particle path (\int f(r) ds or \int f(r) dt )
      !------------------------------------------------------------------
      IF(.NOT. ReverseTime .OR. Particles % DtSign == 1 ) THEN
        CALL ParticlePathIntegral( Particles, istep )
      END IF
        
      InitTimestep = .FALSE.
    END DO 

    NoMoving = Particles % NumberOfMovingParticles
    NoMoving = ParallelReduction( NoMoving ) 
    WRITE (Message,'(A,I0,A,I0,A)') 'Timestep ',i,' with ',NoMoving,' moving particles'
    CALL Info(Caller,Message,Level=6)

    ! Freeze particles that are too old
    !----------------------------------------------------------------
    CALL SetRetiredParticles( )    
    
    IF( InfoActive( 15 ) ) THEN 
      CALL ParticleInformation(Particles, ParticleStepsTaken, &
          TimeStepsTaken, tottime )
    END IF
      
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
  
  IF( ReverseTime ) THEN
    IF( Particles % DtSign == -1 ) THEN
      CALL Info(Caller,'Reversing time and restarting integration', Level=7)
      Particles % DtSign = 1
      InitTimestep = .TRUE.
      IF( ASSOCIATED( PTimeVar ) ) THEN
        PTimeVar % Values = 0.0_dp
      END IF      
      DO i = 1, Particles % NumberOfParticles
        CALL SetParticleStatus( Particles, i, PARTICLE_MOVING)
      END DO
      GOTO 1
    ELSE
      CALL Info(Caller,'Time reversal finished!',Level=7)
    END IF
  END IF
  
  IF( ListGetLogical( Params,'Fix Tangent Velocity',GotIt ) ) THEN
    CALL ApplyTangentFix(.TRUE.)
  END IF

  
  CALL Info(Caller,'All done',Level=4)
  CALL Info(Caller, '-----------------------------------------', Level=4 )
  

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
    INTEGER :: FixedCount
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
    FixedCount = 0
    
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
          IF( Particles % Status( k ) == PARTICLE_FIXEDCOORD ) CYCLE

          Particles % Status( k ) = PARTICLE_FIXEDCOORD
          FixedCount = FixedCount + 1
	  IF(.NOT. GotTime ) CYCLE

	  IF( Found ) THEN
            PTimeVar % Values( k ) = PTime(i)
          ELSE
            PTimeVar % Values( k ) = PTimeConst
          END IF 
        END DO
      END DO
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
          IF( Particles % Status( k ) == PARTICLE_FIXEDCOORD ) CYCLE

          Particles % Status( k ) = PARTICLE_FIXEDCOORD
          FixedCount = FixedCount + 1
          IF(.NOT. GotTime ) CYCLE

	  IF( Found ) THEN
            PTimeVar % Values( k ) = PTime(i)
          ELSE
            PTimeVar % Values( k ) = PTimeConst
          END IF 
        END DO
      END DO
    END IF
    DEALLOCATE( PTime, PCond )
              
    CALL Info('SetFixedParticles','Number of new fixed particles: '&
        //I2S(FixedCount),Level=7)
    
  END SUBROUTINE SetFixedParticles 

  
  !> Eliminate particles that are too old sit an fixed boundaries.
  !-------------------------------------------------------------------------
  SUBROUTINE SetRetiredParticles( )

    INTEGER :: i,j,k
    REAL(KIND=dp) :: MaxIntegTime
    LOGICAL :: Found, GotMaxIntegTime
    INTEGER :: FixedCount
    TYPE(ValueList_t), POINTER :: Params

    Params => GetSolverParams()

    MaxIntegTime = ListGetCReal( Params,'Max Integration Time',GotMaxIntegTime )
        
    IF( .NOT. GotMaxIntegTime ) RETURN		  

    PtimeVar => ParticleVariableGet( Particles, 'particle time' )
    IF(.NOT. ASSOCIATED( PTimeVar ) ) THEN
      CALL Fatal('SetRetiredParticles','Cannot retire particles without time!')
    END IF

    FixedCount = 0    
    DO k=1,SIZE( PTimeVar % Values ) 
      IF( Particles % Status( k ) == PARTICLE_FIXEDCOORD ) CYCLE
      IF( PTimeVar % Values( k ) >= MaxIntegTime ) THEN
        FixedCount = FixedCount + 1
        Particles % Status( k ) = PARTICLE_FIXEDCOORD
      END IF
    END DO

    IF( FixedCount > 0 ) THEN
      CALL Info('SetRetiredParticles','Number of new retired particles: '&
          //I2S(FixedCount),Level=7)    
    END IF
      
  END SUBROUTINE SetRetiredParticles


  SUBROUTINE ApplyTangentFix(UndoFix)
    LOGICAL :: UndoFix
    INTEGER :: sgn,i,jv,jf,vdofs
    TYPE(Variable_t), POINTER :: VeloVar, FixVar
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName

    sgn = -1
    IF(UndoFix) sgn=1
    
    VariableName = ListGetString(Params,'Velocity Variable Name',Found)
    IF(.NOT. Found) VariableName = 'flow solution'
    VeloVar => VariableGet( Mesh % Variables, TRIM(VariableName) )    

    FixVar => VariableGet( Mesh % Variables,'fixvelo')
    IF(.NOT. ASSOCIATED(FixVar)) CALL Fatal('ApplyTangentFix','"FixVelo" not found!')
    
    vdofs = VeloVar % dofs
    
    DO i=1, Mesh % NumberOfNodes
      jf = FixVar % Perm(i)
      jv = VeloVar % Perm(i)
      IF(jf==0 .OR. jv==0) CYCLE
      
      VeloVar % Values(vdofs*(jv-1)+3) = VeloVar % Values(vdofs*(jv-1)+3) +&
          sgn * FixVar % Values(jf)
    END DO
    
  END SUBROUTINE ApplyTangentFix


  
  !------------------------------------------------------------------------
  !> Compute field values at the given points in the FE mesh. 
  !-------------------------------------------------------------------------
  SUBROUTINE SetParticleVelocities( FirstStep )
    LOGICAL :: FirstStep
    
    TYPE(Element_t), POINTER :: BulkElement
    INTEGER :: No, Status
    REAL(KIND=dp) :: Coord(3),Velo(3),GradVelo(3,3)    
    TYPE(Element_t), POINTER :: BulkElement2
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Valuelist_t), POINTER :: Params
    REAL(KIND=dp) :: VeloAtPoint(3), GradVeloAtPoint(3,3),dtime
    LOGICAL :: Stat, UseGradVelo, Visited = .FALSE.
    INTEGER :: i,j,k,l,n,TimeOrder, NewLost(3), OldLost, FixedLost
    INTEGER, POINTER :: NodeIndexes(:), FieldPerm(:),FieldIndexes(:)
    REAL(KIND=dp) :: SqrtElementMetric, Weight, Speed, SpeedMin
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:), Coordinate(:,:), Velocity(:,:)
    LOGICAL :: GotIt, SkipZeroTime
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName
    TYPE(Variable_t), POINTER :: VeloVar
    TYPE(Variable_t), POINTER :: DtVar	
    
    
    SAVE :: Visited, Mesh, Basis, dBasisdx, Params, VeloVar, UseGradvelo, DtVar, &
	SpeedMin, NewLost 

    IF( .NOT. Visited ) THEN
      Mesh => GetMesh()
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

    SkipZeroTime = .NOT. ( Particles % DtConstant .OR.  FirstStep ) 
    
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
        Particles % Status(No) = PARTICLE_LOST
        NewLost(1) = NewLost(1) + 1
        CYCLE       
      END IF

      ! If the particle has not moved then it cannot have
      ! any change in the velocity.
      IF( SkipZeroTime ) THEN
        IF( ABS( DtVar % Values(No) ) < TINY( dtime ) ) CYCLE
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
    INTEGER :: i,j,k,l,n,nsize,wallcount,NoVar,NoNorm,dofs,maxdim,VarType
    INTEGER, POINTER :: NodeIndexes(:), PPerm(:)
    REAL(KIND=dp) :: SqrtElementMetric, Norm, PrevNorm = 0.0_dp, Change
    REAL(KIND=dp), POINTER :: Basis(:)
    LOGICAL :: GotIt, Difference,Cumulative,Derivative,GotVar,GotRes,GotOper,&
        UsePerm,InternalVariable,Initiated, Parallel, InitIntegral
    CHARACTER(LEN=MAX_NAME_LEN) :: VariableName, ResultName, OperName, Name
    TYPE(Variable_t), POINTER :: TargetVar, ResultVar, DataVar, Var
    TYPE(Variable_t), POINTER :: ParticleVar
    REAL(KIND=dp), POINTER :: TmpValues(:), NodeValues(:), NewValues(:)
    INTEGER, POINTER :: TmpPerm(:), UnitPerm(:)
    REAL(KIND=dp) :: x,y,z
    
    
    SAVE :: Visited, PrevNorm, UnitPerm

    IF( ReverseTime ) THEN
      IF( Particles % DtSign == -1 ) THEN
        CALL Info(Caller,'Setting the advected fields at turning point',Level=10)
      ELSE
        CALL Info(Caller,'Setting the advected fields after reverse tracking',Level=10)
      END IF
    ELSE
      CALL Info(Caller,'Setting the advected fields',Level=10)
    END IF
    
    Mesh => GetMesh()
    maxdim = dim
    n = Mesh % MaxElementNodes
    ALLOCATE( Basis(n) )
    Coord = 0.0_dp
    Velo = 0.0_dp
    
    Params => GetSolverParams()
    NoNorm = GetInteger( Params,'Norm Variable Index',GotIt)
    NoParticles = Particles % NumberOfParticles
    
    DataVar => VariableGet( Mesh % Variables,'AdvectorData')
    IF( ASSOCIATED( DataVar ) ) THEN
      nsize = SIZE( DataVar % Values )
      DataVar % Output = .FALSE.
    ELSE
      nsize = Mesh % NumberOfNodes
    END IF    

    Parallel = ( ParEnv % PEs > 1 ) 

    Initiated = .FALSE.
100 NoVar = 0 
    DO WHILE(.TRUE.)
      NoVar = NoVar + 1
      
      WRITE (Name,'(A,I0)') 'Variable ',NoVar
      VariableName = GetString( Params,Name,GotVar)
      IF(.NOT. GotVar ) EXIT
      
      WRITE (Name,'(A,I0)') 'Operator ',NoVar
      OperName = GetString( Params,Name,GotOper)

      InitIntegral = .FALSE.
      IF( ReverseTime ) THEN
        IF( VariableName == 'particle disp' .OR. VariableName == 'particle disp abs' ) THEN
          ! This is done only at end of forward cycle
          IF( Particles % DtSign == -1 ) CYCLE          
        ELSE IF( VariableName == 'particle time integral' .OR. &
            VariableName == 'particle distance integral' ) THEN
          ! These are done both at forward and backward cycle
          ! End of backward cycle, initialize the integral
          ! End of forward cycle, project particle info to nodal info
          InitIntegral = ( Particles % DtSign == -1 )
        ELSE
          ! This is done at end of backward cycle
          IF( Particles % DtSign == 1 ) CYCLE
        END IF
      END IF
      
      CALL Info(Caller,'Setting field '//I2S(NoVar)//' for variable: '//TRIM(VariableName),Level=15)

      ! Get the target variables
      ! Variables starting with 'particle' as associated with particles
      !----------------------------------------------------------------
      TargetVar => NULL()
      ResultVar => NULL()
      
      IF( VariableName == 'particle coordinate' .OR. &
          VariableName == 'particle velocity' .OR. &
          VariableName == 'particle force' .OR. &
          VariableName == 'particle disp') THEN

        dofs = dim 
        InternalVariable = .TRUE.
      ELSE IF( SEQL(VariableName, 'particle') ) THEN
        dofs = 1
        InternalVariable = .TRUE.
      ELSE
        TargetVar => VariableGet( Mesh % Variables, TRIM(VariableName) )
        IF( .NOT. ASSOCIATED(TargetVar)) THEN
          CALL Fatal(Caller,'Could not obtain target variable:'&
              //TRIM(VariableName))
          CYCLE
        END IF
        dofs = TargetVar % dofs
        IF( dofs > 1 ) THEN
          CALL Fatal(Caller,'Advection implemented so far only for scalars')
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
      
      CALL Info(Caller,'Using operator '//I2S(NoVar)//' for variable: '//TRIM(OperName),Level=15)
      
      WRITE (Name,'(A,I0)') 'Result Variable ',NoVar
      ResultName = GetString( Params,Name,GotRes)
      IF( .NOT. GotRes ) THEN
        ResultName = TRIM(OperName)//' '//TRIM(VariableName)
      END IF


      ! Create variables if they do not exist
      !---------------------------------------------------------      
      ResultVar => VariableGet( Mesh % Variables, TRIM(ResultName) )
      IF( ASSOCIATED(ResultVar) ) THEN        
        IF( ASSOCIATED( DataVar) ) THEN        
          IF( DataVar % TYPE /= ResultVar % TYPE ) THEN
            CALL Fatal(Caller,'ResultVar is of wrong type, use new name for result variable!')
          END IF
          IF( SIZE( DataVar % Values ) /= SIZE( ResultVar % Values) ) THEN
            CALL Fatal(Caller,'ResultVar is of wrong size, use new name for result variable!')
          END IF
        END IF
        CALL Info(Caller,'Found a pre-existing result variable: '//TRIM(ResultName),Level=20)
      ELSE
        GotIt = .FALSE.
        PPerm => NULL()
        IF( ASSOCIATED( DataVar ) ) THEN
          CALL Info(Caller,'Using non-nodal given permutation for data',Level=15)
          PPerm => DataVar % Perm
          VarType = DataVar % TYPE
        ELSE IF( ASSOCIATED( TargetVar ) ) THEN
          CALL Info(Caller,'Using inherited permutation for data',Level=15)
          PPerm => TargetVar % Perm
          VarType = TargetVar % TYPE
        END IF

        IF( .NOT. ASSOCIATED( PPerm ) ) THEN
          IF(.NOT. ASSOCIATED(UnitPerm)) THEN
            CALL Info(Caller,'Creating unity permutation for data',Level=15)
            ALLOCATE( UnitPerm(nsize) )
            DO i=1,nsize
              UnitPerm(i) = i
            END DO
          END IF
          PPerm => UnitPerm
          VarType = 0
        END IF
        
        CALL VariableAddVector( Mesh % Variables,Mesh,PSolver,ResultName,dofs,&
            Perm = PPerm, VarType = VarType )
        
        IF( dofs == 1 ) THEN
          CALL Info(Caller,'Created a scalar variable: '//TRIM(ResultName) )
        ELSE
          CALL Info(Caller,'Created a vector variable: '//TRIM(ResultName) )          
        END IF
        ResultVar => VariableGet( Mesh % Variables, TRIM(ResultName))
        IF(.NOT. ASSOCIATED(ResultVar)) CALL Fatal(Caller,'Problems in VariableAdd')
      END IF


      ! Finally, set the values
      !---------------------------------------------------------      
      IF( InternalVariable .AND. .NOT. InitIntegral ) THEN
        CALL Info(Caller,'Setting particle variable "'//TRIM(VariableName)//'" to fields',Level=15)

        IF( VariableName == 'particle disp abs') THEN
          IF( ResultVar % Dofs /= 1 ) THEN
            CALL Fatal(Caller,'Variable should have one dof: '//TRIM(VariableName))
          END IF
          DO i=1,NoParticles
            NewValues(i) = SQRT(SUM((Particles % Coordinate(i,:) - OrigCoordinate(i,:))**2))
          END DO

        ELSE IF( VariableName == 'particle disp') THEN
          IF( ResultVar % Dofs /= dim ) THEN
            CALL Fatal(Caller,'Variable should have dim dofs: '//TRIM(VariableName))
          END IF
          DO i=1,NoParticles
            DO j=1,dim
              NewValues(dim*(i-1)+j) = Particles % Coordinate(i,j) - OrigCoordinate(i,j)
            END DO
          END DO

        ELSE IF( VariableName == 'particle coordinate') THEN
          IF( ResultVar % Dofs /= dim ) THEN
            CALL Fatal(Caller,'Variable should have dim dofs: '//TRIM(VariableName))
          END IF          
          DO j=1,dim
            DO i=1,NoParticles
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
            CALL Fatal(Caller,'Variable should have dim dofs: '//TRIM(VariableName))
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
            CALL Fatal(Caller,'Variable should have dim dofs: '//TRIM(VariableName))
          END IF
          DO i=1,NoParticles
            DO j=1,dim
              NewValues(dim*(i-1)+j) = Particles % Force(i,j)
            END DO
          END DO

        ELSE IF( VariableName == 'particle status') THEN
          NewValues(1:NoParticles) = 1.0_dp * Particles % Status(1:NoParticles)

        ELSE IF( VariableName == 'particle partition') THEN
          IF( ASSOCIATED( Particles % Partition ) ) THEN
            NewValues(1:NoParticles) = ParEnv % MyPe + 1.0_dp 
          END IF
            
        ELSE IF( VariableName == 'particle number') THEN
          DO i=1,NoParticles
            NewValues(i) = 1.0_dp * i
          END DO

        ELSE IF( VariableName == 'particle index') THEN
          NewValues(1:NoParticles) = 1.0_dp * Particles % NodeIndex(1:NoParticles)
          
        ELSE IF( SEQL(VariableName, 'particle') ) THEN
          ParticleVar => ParticleVariableGet( Particles, VariableName )
          IF( ASSOCIATED( ParticleVar ) ) THEN
            NewValues(1:NoParticles) = ParticleVar % Values(1:NoParticles)
          ELSE
            CALL Warn(Caller,'Field does not exist: '//TRIM(VariableName))
          END IF
        END IF
        
      ELSE 
        CALL Info(Caller,'Setting field variable "'//TRIM(VariableName)//'" to advected fields',Level=15)

        IF( InitIntegral ) THEN
          ! This is a cludge as the path integrals are somewhat special
          TargetVar => ResultVar
        END IF
        
        
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
        
        ! These are initializations related to path integrals.
        ! No parallel communication required for these.      
        IF( InitIntegral ) THEN
          IF( ASSOCIATED( TargetVar ) ) THEN
            ParticleVar => ParticleVariableGet( Particles, VariableName )
            IF( ASSOCIATED( ParticleVar ) ) THEN
              CALL Info(Caller,'Initialize "'//TRIM(ParticleVar % Name)//'" with "'//TRIM(TargetVar % Name)//'"')
              ParticleVar % Values = NewValues(1:NoParticles)
            END IF
          END IF          
          TargetVar => NULL()
          CYCLE
        END IF

      END IF
        
      ! In a serial case the nodes and particles are directly associated. 
      ! In a parallel case we need to transfer the values from particles in 
      ! different partitions to nodes. 
      !---------------------------------------------------------------------
      IF( Parallel ) THEN
        NodeValues = 0._dp
        CALL ParticleAdvectParallel( Particles, NewValues, NodeValues, dofs )
      END IF

      ! Finally move the nodal values to the target variable 
      !---------------------------------------------------------------------
      IF( ASSOCIATED( DataVar ) ) THEN
        n = SIZE(ResultVar % Values)
        IF( Difference .OR. Derivative ) THEN
          ResultVar % Values = NodeValues(1:n) - TargetVar % Values 
        ELSE IF( Cumulative ) THEN
          ResultVar % Values = NodeValues(1:n) + ResultVar % Values
        ELSE
          ResultVar % Values = NodeValues(1:n)
        END IF
      ELSE
        DO i=1,nsize
          k = i
          IF( ASSOCIATED( ResultVar % Perm ) ) k = ResultVar % Perm( k )
          IF( k == 0 ) CYCLE
          dofs = ResultVar % Dofs
          DO j=1,dofs
            IF( Difference .OR. Derivative ) THEN
              ResultVar % Values( dofs*(k-1)+j ) = NodeValues( dofs*(i-1)+j ) - TargetVar % Values( dofs*(k-1)+j ) 
            ELSE IF( Cumulative ) THEN
              ResultVar % Values( dofs*(k-1)+j ) = NodeValues( dofs*(i-1)+j ) + ResultVar % Values( dofs*(k-1)+j )
            ELSE
              ResultVar % Values( dofs*(k-1)+j ) = NodeValues( dofs*(i-1)+j ) 
            END IF
          END DO
        END DO
      END IF
      
      IF( Derivative ) ResultVar % Values = ResultVar % Values / dertime 

      BLOCK 
        INTEGER :: t, LocalPerm(10)
        REAL(KIND=DP) :: cval
        TYPE(Element_t), POINTER :: Element
        REAL(KIND=dp) :: DgScale
        LOGICAL :: GotScale

        IF( ResultVar % TYPE == variable_on_nodes_on_elements ) THEN
       
          DGScale = ListGetCReal( Solver % Values,'DG Nodes Scale',GotScale )
          IF(.NOT. GotScale ) DgScale = 1.0 / SQRT( 3.0_dp ) 
          GotScale = ( ABS( DGScale - 1.0_dp ) > TINY( DgScale ) )
          
          IF( GotScale ) THEN
            CALL Info(Caller,'Expanding shrinked DG field',Level=12)        
            DO t=1, Mesh % NumberOfBulkElements
              Element => Mesh % Elements(t)
              n = Element % TYPE % NumberOfNodes
              LocalPerm(1:n) = ResultVar % Perm( Element % DGIndexes )
              IF( ANY( LocalPerm(1:n) == 0) ) CYCLE
              cval = SUM( ResultVar % Values( LocalPerm(1:n) ) ) / n
              DO i=1,n
                j = LocalPerm(i)
                ResultVar % Values(j) = cval + ( ResultVar % Values(j) - cval ) * ( 1.0_dp / DgScale )
              END DO
            END DO
          END IF
        END IF
      END BLOCK
      
      ! To allow computation of change in the standard manner the Variable
      ! is set to point to the one of interest. This is mainly used in the 
      ! tests, or could be used in for convergence monitoring also. 
      !---------------------------------------------------------------
      IF( NoVar == NoNorm ) THEN
        n = SIZE( ResultVar % Values ) 
        Norm = SUM( ResultVar % Values ** 2)

        n = ParallelReduction(n)
        Norm = ParallelReduction(Norm)
        Norm = SQRT(Norm/n)
        
        Change = 2.0 * ABS( Norm-PrevNorm ) / ( Norm + PrevNorm )
        PrevNorm = Norm
               
        Solver % Variable % Norm = Norm
        Solver % Variable % NonlinChange = Change
        Solver % Variable % Values = Norm

        ! Here the name is ComputeChange in order to get the change also to ElmerGUI
        ! albeit in a very dirt style. One could also edit ElmerGUI....
        WRITE( Message, '(a,g15.8,g15.8,a)') &
            'SS (ITER=1) (NRM,RELC): (',Norm, Change,' ) :: '//TRIM( ResultName )
        CALL Info( 'ComputeChange', Message, Level=3 )
      END IF
    END DO

    ! Allocate the local new temporal values
    IF(.NOT. Initiated ) THEN
      CALL Info(Caller,'Allocating for temporal value vectors',Level=15)
      NoVar = NoVar - 1
      IF( NoVar < 1 ) THEN
        CALL Fatal(Caller,'No target and result variables exist!')
      END IF
      ALLOCATE( NewValues( maxdim * Particles % NumberOfParticles ) ) 
      NewValues = 0.0_dp
      IF( Parallel ) THEN
        ALLOCATE( NodeValues( maxdim * nsize ) )
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
  LOGICAL :: Found, AdvectElemental, AdvectDG, AdvectIp
  INTEGER :: NormInd

  Params => GetSolverParams()

  ! These are default setting that make the operation of the advection solver 
  ! possible. There should always be one passive particle for each active node.
  !---------------------------------------------------------------------------
  AdvectElemental = ListGetLogical( Params,'Advect Elemental',Found) 
  AdvectDG = ListGetLogical( Params,'Advect DG',Found) 
  AdvectIp = ListGetLogical( Params,'Advect Ip',Found) 
  
  IF( AdvectElemental ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-elem AdvectorData')
  ELSE IF( AdvectDg ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-dg AdvectorData')
  ELSE IF( AdvectIp ) THEN
    CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-ip AdvectorData')
    ! This is not a good idea since if the solver is not active this will have size zero!
    !ELSE
    !  CALL ListAddString( Params,NextFreeKeyword('Exported Variable',Params),'-nodal AdvectorData')
  END IF
    
  CALL ListAddInteger( Params,'Time Order',0 )
  CALL ListAddNewLogical( Params,'Particle Accurate At Face',.FALSE.)  

  ! If we want to show a pseudonorm add a variable for which the norm
  ! is associated with.
  CALL ListAddNewString( Solver % Values,'Variable','-nooutput -global particleadvector_var')
  CALL ListAddNewLogical( Params,'No Matrix',.TRUE.)
  
END SUBROUTINE ParticleAdvector_Init
