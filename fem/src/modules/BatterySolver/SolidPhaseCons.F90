!------------------------------------------------------------------------------
! For copyrights see the BatteryUtils.F90 file in this directory!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Initialization of the primary solver. 
!> If requested create an internal mesh.
!> This is the only solver that operates on the internal 1D mesh
!> that is active in every node of the global 1D/2D/3D mesh. 
!------------------------------------------------------------------------------
SUBROUTINE SolidPhaseCons_Init( Model,Solver,dt,Transient)
  USE DefUtils
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient, Found
  !------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'SolidPhaseCons_init'
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Mesh_t), POINTER :: Mesh, PMesh

  Params => GetSolverParams()

  CALL ListAddInteger(Params,'1D Active Direction',1)
  CALL ListAddNewConstReal(Params,'1D Mesh Length',1.0_dp)
  
  ! Create 1D mesh on the fly and set it to be the active mesh of this solver.
  IF( GetLogical( Params,'1D Mesh Create') ) THEN
    CALL Info(Caller,'Creating internal 1D mesh')

    ! The initial 1D mesh is x \in [0,1].
    ! Currently we deal with the unite mesh inside the code
    Mesh => CreateLineMesh( Params )
    
    Mesh % OutputActive = .FALSE.
    Solver % Mesh => Mesh 

    
    ! Add the mesh to the list of meshes
    PMesh => Model % Meshes
    IF( ASSOCIATED( PMesh ) ) THEN
      DO WHILE ( ASSOCIATED( PMesh % Next ) ) 
        Pmesh => PMesh % Next
      END DO
      Pmesh % Next => Mesh
    END IF
  END IF

  ! This variable is used to solve one stride at a time.
  CALL ListAddNewString( Params,'Variable','Cs OneDim')

  ! This is 1D mesh: never optimize the bandwidth so we know it is identity!
  CALL ListAddLogical( Params,'Optimize Bandwidth',.FALSE.)

  ! We will assembly mass & stiffness matrices only once and use them
  ! for timestepping. To be able to separate mass matrix we need to have this flag on. 
  CALL ListAddLogical( Params,'Use Global Mass Matrix',.TRUE.)

  ! We solve number of 1D equation but study the average norm rather than an individual one
  CALL ListAddNewLogical( Params,'Skip Compute Nonlinear Change',.TRUE.)
  CALL ListAddNewLogical( Params,'Skip Compute Steady State Change',.TRUE.)
  
  IF( ListGetLogical( Params,'Linearize Flux',Found ) ) THEN
    CALL ListAddNewLogical( Params,'Calculate Cs Sensitivity',.TRUE.)
  END IF
  
END SUBROUTINE SolidPhaseCons_Init



!------------------------------------------------------------------------------
!> Solve for the 1D diffusion equation taking place at every node (or just once).
!> By default the equation is solved in a spherically symmetric domain (dim=3).
!> Also line strip (dim=1) and circle (dim=2) could be possible geometries but
!> are not currently supported. 
!------------------------------------------------------------------------------
SUBROUTINE SolidPhaseCons( Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  USE DefUtils
  USE BatteryModule
  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t):: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  CHARACTER(*), PARAMETER :: Caller = 'SolidPhaseCons'
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Matrix_t),POINTER  :: Amat
  TYPE(Element_t),POINTER :: Element, Element1D
  TYPE(ValueList_t), POINTER :: Material, PrevMaterial
  TYPE(Variable_t), POINTER :: Var, SensVar
  INTEGER :: t,i,j,k,l,n,nn,dofs,elems,iLeft,iRight,jLeft,jRight,&
      CsNodes, timestep, prevtimestep = -1, iter, maxiter, &
      VisitedTimes = 0, NoLimited, NoPassive
  INTEGER, POINTER :: CsPerm(:)  
  LOGICAL, ALLOCATABLE :: NodeDone(:)
  LOGICAL :: Found, Newton, Show, DoRelax, DoSSRelax, &
      LimitDxLoc, LimitDxGlo, NewtonConst, DoPotRelax
  REAL(KIND=dp), POINTER :: x(:), xprev(:,:), bvec(:)
  REAL(KIND=dp), ALLOCATABLE :: x0(:), Fullx0(:), phis0(:), phie0(:)
  REAL(KIND=dp) :: DiffCoeff, FluxCoeff, Norm, SumNorm, dNorm, SumdNorm, &
      Change, NonlinTol, Flux, DiffFlux, R_s, Relax, SSRelax, TotFlux0(2), &
      dx, dxmax, dxlimg, dxliml, NewtonFactor, NewtonCoeff, SSRelax0, &
      PrevSSRelax, PotRelax
  CHARACTER(LEN=MAX_NAME_LEN) :: str
  LOGICAL :: Visited = .FALSE., DoOutput, PostFix
  REAL(KIND=dp) :: fluxerr, minfluxerr, maxfluxerr
  
  SAVE Visited, prevtimestep, PrevSSRelax, &
      CsNodes, CsPerm, iLeft, iRight, jLeft, jRight, &
      x0, Fullx0, Phie0, Phis0

  !------------------------------------------------------------------------------

  CALL Info(Caller,'-------------------------------------------------------')
  CALL Info(Caller,'Solving 1D evolution equation for each solid phase node')
  CALL Info(Caller,'-------------------------------------------------------')

  CALL InitializeBattery()


  Params => GetSolverParams()

  SubMesh => GetMesh()
  dofs = SubMesh % NumberOfNodes
  elems = SubMesh % NumberOfBulkElements

  ! Variable associated to 1D mesh
  Var => Solver % Variable

  ! Set pointers to the already allocated matrix equation: Ax=b
  Amat => Solver % Matrix
  bvec => Amat % Rhs
  x => Var % Values
  xprev => Var % PrevValues

  maxiter = ListGetInteger( Params,'Nonlinear System Max Iterations',Found)
  IF(.NOT. Found ) maxiter = 1  
  Relax = ListGetCReal( Params,'Nonlinear System Relaxation Factor',DoRelax )

  SSRelax0 = ListGetCReal( Params,'Solid Phase Relaxation Factor',DoSSRelax )
  IF(.NOT. DOSSRelax ) SSRelax0 = 1.0
  PrevSSRelax = SSRelax0
  
  NonlinTol = ListGetCReal( Params,'Nonlinear System Convergence Tolerance',Found )

  dxliml = dt * ListGetCReal( Params,'Maximum Local Change Speed',LimitDxLoc )  
  dxlimg = dt * ListGetCReal( Params,'Maximum Global Change Speed',LimitDxGlo )  
  
  
  IF( .NOT. Visited ) THEN
    CALL Info(Caller,'Local mesh has '//I2S(dofs)//' nodes and '//I2S(elems)//' elements')

    IF( SIZE( Var % Values ) /= dofs ) THEN
      CALL Fatal(Caller,'Size of variable "'//TRIM(Var % Name)//&
          '" should be equal to number of 1D mesh nodes!')
    END IF
    
    ! We assume simple mesh that extends in the 1st coordinate direction.
    ! The BCs can only be given at the start and at the finish of it. 
    iLeft = ExtremeLeftNode( SubMesh ) 
    iRight = ExtremeRightNode( SubMesh ) 
    
    jLeft = Var % Perm(iLeft)
    jRight = Var % Perm(iRight)    
    
    ! The 1D solid phase diffusion equation is solved in SubMesh.
    ! The other equations are solved in MainMesh.
    !-------------------------------------------------------------------------------------
    CsPerm => CsVar % Perm  
    ! Number Of active nodes for concentration
    CsNodes = SIZE( CsVar % Values )    
    
    CALL Info(Caller,'Allocating full solution for solid phase',Level=8)

    ! If requested saves the solid phase profile for visualization
    ! This cannot be created before since we dont easily know size of 1D mesh.
    DoOutput = ListGetLogical( Params,'Save Solid Phase Profile', Found )
    CALL VariableAddVector( MainMesh % Variables, MainMesh, Solver,'Cs Profile',dofs,&
        Perm = CsVar % Perm, Output = DoOutput, Secondary = .TRUE.)
    
    CsFullVar => VariableGet( MainMesh % Variables,'Cs Profile' )
    n = SIZE( Var % PrevValues, 2 )

    ALLOCATE( CsFullVar % PrevValues(dofs*CsNodes, n ) )
    ! Copy initial conditions
    DO i=1,SIZE(CsVar % Values)
      CsFullVar % Values((i-1)*dofs+1:i*dofs) = CsVar % Values(i)
    END DO
    DO i=1,n
      CsFullVar % PrevValues(:,i) = CsFullVar % Values(:)
    END DO
    
    ALLOCATE( x0( SIZE( x ) ) )
    IF( DoSSRelax .OR. LimitDxGlo ) ALLOCATE( FullX0( SIZE( CsFullVar % Values ) ) )
    !PrevSSRelax = SSRelax0

    ! Memorize the initial solid phase consentration
    IF( ASSOCIATED( CsInitVar ) ) THEN
      CsInitVar % Values = CsVar % Values
    END IF
    
    Visited = .TRUE.
  END IF
    
  VisitedTimes = VisitedTimes + 1

  NoPassive = GetInteger( Params,'Number of Passive Visits',Found )

  PotRelax = ListGetCReal( Params,'Potential Relaxation Factor',DoPotRelax)
  IF( DoPotRelax ) THEN    
    IF( .NOT. ALLOCATED( PhiS0 ) ) THEN
      ALLOCATE( PhiS0( SIZE( PhisVar % Values ) ), Phie0( SIZE( PhieVar % Values ) ) )
    END IF
    Phis0 = PhisVar % Values
    Phie0 = PhieVar % Values
    PhisVar % Values = PotRelax * PhisVar % Values + (1.0_dp-PotRelax) * Phis0
    PhieVar % Values = PotRelax * PhieVar % Values + (1.0_dp-PotRelax) * Phie0    
  END IF
      
  IF( VisitedTimes < NoPassive ) THEN
    RETURN
  END IF

  PostFix = GetLogical( Params,'Check Material Balance',Found )
  minfluxerr = HUGE(minfluxerr)
  maxfluxerr = -HUGE(maxfluxerr)
  
  timestep = GetTimestep()

  ! If we are on a new timestep advance concentration forward in time.
  IF( timestep /= prevtimestep ) THEN  
    IF( timestep > 1 ) THEN
      CALL Info(Caller,'Moving concentration forward in time: '//I2S(timestep))
      n = SIZE( CsFullVar % PrevValues,2 )
      DO i=n,2,-1
        CsFullVar % PrevValues(:,i) = CsFullVar % PrevValues(:,i-1)
      END DO
      CsFullVar % PrevValues(:,1) = CsFullVar % Values
    END IF
    prevtimestep = timestep
  END IF

  ALLOCATE( NodeDone( MainMesh % NumberOfNodes ) )
  
  ! Update surface flux
  Newton = ListGetLogical( Params,'Linearize Flux',Found )     
  NewtonConst = ListGetLogical( Params,'Linearize Flux Average',Found)
  NewtonCoeff = ListGetCReal( Params,'Linearize Flux Multiplier',Found)
  IF(.NOT. Found ) NewtonCoeff = 1.0
  IF( NewtonConst ) Newton = .TRUE.
    
  IF( DoSSRelax .OR. LimitDxGlo ) Fullx0 = CsFullVar % Values 

#if 0
  TotFlux0(1) = SUM( JliVar % Values * AnodeWeight, AnodeWeight > 0 )
  TotFlux0(2) = SUM( -JliVar % Values * AnodeWeight, AnodeWeight < 0 )
#endif
  
  DO iter = 1, maxiter
    IF(maxiter>1) CALL Info(Caller,'Nonlinear system iteration: '//I2S(iter),Level=5)
    NoLimited = 0

#if 0
    IF( iter > 1 .AND. ListGetLogical( Params,'Fix Potential', Found ) )  THEN
      BLOCK
        REAL(KIND=dp) :: TotSens(2), TotFlux(2), dphi(2)
        INTEGER :: iter
        
        CALL ListAddLogical( Params,'Calculate Phis Sensitivity',.TRUE.)

        DO iter = 1, 20
        
          CALL ButlerVolmerUpdate( Solver )
          SensVar => VariableGet( MainMesh % Variables,'dJli dPhis')
          IF(.NOT. ASSOCIATED( SensVar ) ) THEN
            CALL Fatal(Caller,'Could not find variable "dJli dPhis"')
          END IF

          TotFlux(1) = SUM( JliVar % Values * AnodeWeight, AnodeWeight > 0 )
          TotFlux(2) = SUM( -JliVar % Values * AnodeWeight, AnodeWeight < 0 )

          TotSens(1) = SUM( SensVar % Values * AnodeWeight, AnodeWeight > 0 )
          TotSens(2) = SUM( -SensVar % Values * AnodeWeight, AnodeWeight < 0 )

          ! TotFlux + TotSens * dPhi = TotFLux0
          dphi = 0.5 * ( TotFlux0 - TotFlux ) / TotSens

          PRINT *,'TotSens:',iter,TotSens,TotFlux,dPhi
          WHERE( AnodeWeight > 0 )
            PhisVar % Values = PhisVar % Values + dphi(1)
          ELSE WHERE
            PhisVar % Values = PhisVar % Values + dphi(2)
          END WHERE
        END DO
        !CALL ListAddLogical( Params,'Calculate Phis Sensitivity',.FALSE.)        
      END BLOCK
    END IF
#endif
    
    IF( iter == 1 ) THEN    
      CALL ButlerVolmerUpdate( Solver )
      ! On the 1st iteration save the flux 
      IF( UseMeanFlux ) THEN
        Jli0 = JliVar % Values
      END IF
    END IF
      
    IF( Newton ) THEN
      SensVar => VariableGet( MainMesh % Variables,'dJli dCs')
      IF(.NOT. ASSOCIATED( SensVar ) ) THEN
        CALL Fatal(Caller,'Could not find variable "dJli dCs"')
      END IF
      IF( NewtonConst ) THEN
        NewtonFactor = SUM( ABS( SensVar % Values ) ) / SIZE( SensVar % Values )
        PRINT *,'newton factor:',NewtonFactor, NewtonCoeff, &
            MINVAL( SensVar % Values), MAXVAL( SensVar % Values )
      END IF
    END IF

    SumNorm = 0.0_dp
    SumdNorm = 0.0_dp
    PrevMaterial => NULL()
    NodeDone = .FALSE.
    dxmax = 0.0_dp
    
    ! Loop over elements is done only in order to have handle for Material
    DO t=1,MainMesh % NumberOfBulkElements
      Element => MainMesh % Elements(t)
      IF( ANY( CsPerm(Element % NodeIndexes) == 0) ) CYCLE

      Material => GetMaterial(Element)

      n  = GetElementNOFNodes(Element)
      nn = GetElementNOFDOFs(Element)
      
      DO i=1,n
        j = Element % NodeIndexes(i)
        IF( NodeDone(j) ) CYCLE           
        NodeDone(j) = .TRUE.
        k = CsPerm(j)

        ! Copy the 1D concentration values related to a node        
        x(1:dofs) = CsFullVar % Values(dofs*(k-1)+1:dofs*k)
        xprev(1:dofs,:) = CsFullVar % PrevValues(dofs*(k-1)+1:dofs*k,:)
        x0 = x
        
        ! Obtain flux from precomputed Butler-Volmer solution
        IF( UseMeanFlux .OR. UseTimeAveFlux ) THEN
          Flux = 0.5_dp * ( JliVar % Values(k) + Jli0(k) ) 
        ELSE
          Flux = JliVar % Values(k)
        END IF
        
        IF( Newton ) THEN
          IF( NewtonConst ) THEN
            diffflux = NewtonCoeff * NewtonFactor
          ELSE
            diffFlux = NewtonCoeff * SensVar % Values(k)
          END IF
          IF( UseMeanFlux .OR. UseTimeAveFlux ) diffFlux = diffFlux / 2
        END IF
                
        ! The 1D diffusion matrix is the same unless the diffusion coefficient has changed.
        IF( .NOT. ASSOCIATED( Material, PrevMaterial ) ) THEN
          Show = .TRUE.
          FluxCoeff = SolidFluxScaling( Material )        
          DiffCoeff = ListGetCReal( Material,'Solid Phase Diffusion Coefficient')    
          R_s = ListGetCReal( Material,'Particle Radius')

          ! We assumed mesh is [0,1]
          ! This means we are missing R_s**3 from the matrix equation and 1/R_s from dBasisdx
          ! The 1st is scaled away also in the r.h.s. but the second means we must rescale diffusion
          DiffCoeff = DiffCoeff / R_s**2
          MassMult1D = 0.0_dp
          
          CALL LocalAssembly1D(DiffCoeff)
          PrevMaterial => Material
          CALL CopyBulkMatrix( Amat, BulkMass = .TRUE. ) 
        ELSE
          IF( t == 0 ) THEN
            PRINT *,'Restore:'
            PRINT *,'BulkValues',Amat % BulkValues
            PRINT *,'BulkMassValues',Amat % BulkMassValues
            PRINT *,'BulkRhs',Amat % BulkRhs
          END IF
          CALL RestoreBulkMatrix( Amat )
        END IF

        CALL Default1stOrderTimeGlobal(Solver)
        
        ! Source term as the Neumann BC in (3.4) 
        bvec(jRight) = bvec(jRight) - FluxCoeff * Flux 
        IF( Newton ) THEN
          CALL CRS_AddToMatrixElement( Amat,jRight,jRight, FluxCoeff*diffFlux)
          bvec(jRight) = bvec(jRight) + FluxCoeff * diffFlux * x0(jRight)
        END IF

        !CALL Default1stOrderTimeGlobal(Solver)

        Norm = DefaultSolve()

        IF( PostFix ) THEN
          fluxerr = SUM( MassMult1D * ( xprev(:,1) - x ) ) / ( FluxCoeff * Flux ) - 1         
          maxfluxerr = MAX( maxfluxerr, fluxerr )
          minfluxerr = MIN( minfluxerr, fluxerr )
        END IF
               
        IF( DoRelax ) THEN
          x = Relax * x + (1-Relax) * x0
        END IF
        
        ! The concentration change at the solid-liquid interface
        dx = x(jright) - x0(jright)

        ! Limit concentration change in each individual node separately
        IF( LimitDxLoc ) THEN
          IF( ABS( dx ) > dxliml ) THEN          
            
            NoLimited = NoLimited + 1
            IF( NoLimited < 5 ) THEN
              PRINT *,'Limited Cs:',iter,j,k,x(Jright),dx
              PRINT *,'xs:',x(1:10)
              PRINT *,'xe:',x(dofs-9:dofs)
              PRINT *,'x0s:',x0(1:10)
              PRINT *,'x0e:',x0(dofs-9:dofs)
              PRINT *,'dxs:',x(1:10)-x0(1:10)
              PRINT *,'dxe:',x(dofs-9:dofs)-x0(dofs-9:dofs)
            END IF

            ! Solve the same equation again, now with Dirichlet BCs such
            ! that the change is solution is no more than "dxliml"
            CALL ZeroRow( Amat, jRight )
            CALL CRS_AddToMatrixElement( Amat,jRight,jRight,1.0_dp)
            bvec(jRight) = x0(jright) + SIGN( dxliml, dx ) 
            Norm = DefaultSolve()
            dx = x(jright) - x0(jright)
          END IF
        END IF

        ! Remember the maximum change in concentration
        dxmax = MAX( ABS(dx), dxmax )
        
        ! We skipped internal computation 
        Norm = SUM( x**2 ) 
        dNorm = SUM( (x-x0)**2 )

        ! Compute summed norm for all the solid phase nodes
        SumNorm = SumNorm + Norm
        SumdNorm = SumdNorm + dNorm
        
        ! Back copy the computed node-wise results to the full vector
        CsFullVar % Values(dofs*(k-1)+1:dofs*k) = x
        CsVar % Values(k) = x(jRight)
      END DO
      
      Show = .FALSE.
    END DO
    
    IF( LimitDxLoc ) THEN
      CALL Info(Caller,'Limiter applied '//I2S(NoLimited)//&
          ' times on iteration '//I2S(iter))
    END IF
    
    ! For testing purposes add the summed average norm as the target norm 
    Change = SQRT( SumdNorm / SumNorm )

    Norm = SQRT( SumNorm / CsNodes )    
    Solver % Variable % Norm = Norm 

    ! We imitate here a standard norm output to be able able to study the
    ! output as usual (using "grep", for example). 
    str = ListGetString( Params,'Equation')
    WRITE( Message, '(a,g15.8,g15.8,a)') &
        'NS (ITER='//i2s(iter)//') (NRM,RELC): (',Norm, Change,&
        ' ) :: '// TRIM(str)
    CALL Info( Caller, Message )        

    IF( Change < NonlinTol ) EXIT
  END DO

  ! Perform steady-state relaxation and/or global limiting of the change in diffusion.
  IF( DoSSRelax .OR. LimitDxGlo ) THEN
    SSRelax = SSRelax0

    IF( LimitDxGlo ) THEN
      SSRelax = MIN(SSRelax, dxlimg / dxmax )
      ! Don't increase the relaxation too fast. Otherwise we may introduce toggling.
      ! This is pure heuristics.
      IF( SSRelax > PrevSSRelax ) SSRelax = SQRT( SSRelax * PrevSSRelax )
      PrevSSRelax = SSRelax
    END IF
      
    CsFullVar % Values = SSRelax * CsFullVar % Values + (1-SSRelax) * Fullx0

    ! Update the surface concentrations in solid phase after complete solution & relaxation 
    CsVar % Values = CsFullVar % Values(jright::dofs)
  END IF
    
  n = COUNT( CsVar % Values < 0 )
  IF( n > 0 ) THEN
    CALL Fatal(Caller,'We got '//I2S(n)//' negative concentration!')       
  END IF

  ! We want to visualize a given node with the 1D solution then we may
  ! copy the values to the 1D mesh. As the same mesh is used for all variables
  ! we have to choose one...
  i = ListGetInteger(Params,'Visualize Node Index',Found )
  IF( Found ) THEN
    j = CsPerm(i)
    ! Copy the 1D concentration values related to a visualization node
    x = CsFullVar % Values(dofs*(j-1)+1:dofs*j)
  END IF

  IF( PostFix ) THEN
    WRITE( Message,'(A,2ES12.4)') 'Flux error range in 1D solution:',minfluxerr,maxfluxerr
    CALL Info(Caller,Message,Level=5)
  END IF

  IF( LimitDxGlo ) THEN
    CALL ListAddConstReal( Model % Simulation,'res: concentration relax',SSRelax )
    WRITE( Message,'(A,ES12.4)') 'Solid phase relaxation factor: ',SSRelax
    CALL Info(Caller,Message,Level=5)
  END IF     
  
  CALL Info(Caller,'Solid phase concentration computed',Level=8)

  
CONTAINS


  ! Assembly the 1D equation for solid phase concentration.
  ! We have scaled the system.
  !----------------------------------------------------------------------
  SUBROUTINE  LocalAssembly1D( DiffCoeff )
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(KIND=dp) :: DiffCoeff
    !------------------------------------------------------------------------------
    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:),dBasisdx(:,:)
    REAL(KIND=dp), ALLOCATABLE, SAVE :: MASS(:,:), STIFF(:,:), FORCE(:)
    REAL(KIND=dp) :: weight, DetJ, rad
    LOGICAL :: Stat,Found
    INTEGER :: i,j,k,t,e,p,q,dim,m,allocstat
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t), SAVE :: Nodes
    !------------------------------------------------------------------------------
    
    ! Allocate storage if needed
    IF (.NOT. ALLOCATED(Basis)) THEN
      m = SubMesh % MaxElementDofs
      ALLOCATE(Basis(m), dBasisdx(m,3),&
          MASS(m,m), STIFF(m,m), FORCE(m), STAT=allocstat)      
      IF (allocstat /= 0) THEN
        CALL Fatal(Caller,'Local storage allocation failed')
      END IF
    END IF

    CALL InitializeToZero( Amat, bvec )
    
    DO e = 1, SubMesh % NumberOfBulkElements
      Element1D => SubMesh % Elements(e)   
      IP = GaussPoints( Element1D, Element1D % Type % GaussPoints2 )

      CALL GetElementNodes( Nodes, UElement=Element1D )      
      m = GetElementNOFNodes(Element1D)

      MASS  = 0._dp
      STIFF = 0._dp
      FORCE = 0._dp
        
      DO t=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element1D, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), detJ, Basis, dBasisdx )
        Weight = IP % s(t) * DetJ

        ! particles are assumed to be spheres
        rad = SUM( Nodes % x(1:m)*Basis(1:m) )
        Weight = Weight * rad**2 ! 4*PI neglected
          
        DO p=1,m
          DO q=1,m
            STIFF(p,q) = STIFF(p,q) + Weight * DiffCoeff * SUM( dBasisdx(p,:) * dBasisdx(q,:) )
            MASS(p,q) = MASS(p,q) + Weight * Basis(p) * Basis(q) 
          END DO

          ! We use this just temporarily to integrate the MassMult1D
          FORCE(p) = FORCE(p) + Weight * Basis(p)
        END DO
        
      END DO

      ! This is the vector that can quickly recover total amounts on 1D stride
      MassMult1D(Element1D % NodeIndexes(1:m)) = MassMult1D(Element1D % NodeIndexes(1:m)) + &
          Force(1:m)
      Force = 0.0_dp
        
      ! When we use global mass matrix this only updates MassValues
      IF(Transient) CALL Default1stOrderTime(MASS,STIFF,FORCE,UElement=Element1D)    

      ! This updates Values and rhs
      CALL DefaultUpdateEquations(STIFF,FORCE,UElement=Element1D)
    END DO
          
  END SUBROUTINE LocalAssembly1D
   
END SUBROUTINE SolidPhaseCons
