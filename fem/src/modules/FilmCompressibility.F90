!------------------------------------------------------------------------------
SUBROUTINE FilmCompressibility_init( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model  
  REAL(KIND=dp) :: dt     
  LOGICAL :: Transient 
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params
  
  Params => GetSolverParams()
  CALL ListAddNewInteger( Params,'Time Derivative Order',0)

  ! This is the instantanoues suggested valued for AC field
  CALL ListAddNewString( Params,'Variable','AC inst')

  ! This is the initial value or relaxed valued for the AC field
  CALL ListAddString( Params,NextFreeKeyword('Exported Variable ',Params),'AC' )
    
END SUBROUTINE FilmCompressibility_Init


!------------------------------------------------------------------------------
!> Subroutine for computing the artificial compressibility from the volume change 
!> of elements. The volume change is obtained by extending the displacement field
!> of a test load to the fluid domain.
!> \ingroup Solvers
!------------------------------------------------------------------------------
SUBROUTINE FilmCompressibility( Model,Solver,dt,Transient )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver   !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model     !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt        !< Timestep size for time dependent simulations
  LOGICAL :: Transient !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Variable_t), POINTER :: HeightSol, PresSol, SensSol, Var, LSVar
  TYPE(Element_t), POINTER :: Element
  TYPE(Mesh_t), POINTER :: Mesh  
  CHARACTER(LEN=MAX_NAME_LEN) :: VarName
  LOGICAL :: FullMode, SensMode, DiffMode, GotEqPres, SolvePDE, LSInit  
  REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:), ElemHeight(:), ElemPres(:), &
      PrevElemHeight(:), PrevElemPres(:), ElemSens(:), Density(:), Thickness(:), EqPres(:)
  REAL(KIND=dp), POINTER :: ACinst(:), ACave(:), gWork(:,:)
  REAL(KIND=dp) :: Norm, ReferencePressure, Diff, dp0, ac_min, &
      RelaxP, RelaxM, Relax, val, Volume, Area, PresInt, grav, rho, &
      dpmax, dhmax, H3Int, LSMin
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: t, i, j, k, n, m, istat, dim, nlift, nVisited = 0
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Nodes_t) :: Nodes 
  LOGICAL :: AllocationsDone = .FALSE., Found, PressureExists, DoIt
  CHARACTER(*), PARAMETER :: Caller = 'FilmCompressibility'

  
  SAVE STIFF, FORCE, Nodes, ElemPres, ElemHeight, ElemSens, Density, Thickness, EqPres, &
      PrevElemHeight, PrevElemPres, ACave, ACinst, AllocationsDone, nVisited 

  CALL Info(Caller,' ')
  CALL Info(Caller,'----------------------------------------------')
  CALL Info(Caller,'Solving compressibility field for FilmPressure')
  CALL Info(Caller,'----------------------------------------------')

  
!------------------------------------------------------------------------------
! Get variables needed for solution
!------------------------------------------------------------------------------
  CALL DefaultStart()
  
  Mesh => Solver % Mesh  
  IF(.NOT. Transient) THEN
    CALL Fatal(Caller,'Implemented only for transient problems!')
  END IF  
  dim = CoordinateSystemDimension()

  ! Instatations suggestion for the AC field value.
  ACinst => Solver % Variable % Values
  
  Params => GetSolverParams()
  
  DiffMode = ListGetLogical( Params,'Differential Mode',Found )
  SensMode =  ListGetLogical( Params,'Sensitivity Mode',Found )
  FullMode = .NOT. (DiffMode .OR. SensMode)
  
  VarName = GetString( Params,'Pressure Variable Name', Found )
  IF ( .NOT. Found ) VarName = 'Pressure'
  PresSol => VariableGet( Mesh % Variables, VarName )
  IF(.NOT. ASSOCIATED(PresSol)) THEN
    CALL Fatal(Caller,'"Pressure" variable '//TRIM(VarName)//' not found!')
  END IF

  VarName = GetString( Params,'Height Variable Name', Found )
  IF ( .NOT. Found ) VarName = 'Height'
  HeightSol => VariableGet( Mesh % Variables, VarName )
  IF(.NOT. ASSOCIATED(HeightSol)) THEN
    CALL Fatal(Caller,'"Height" variable '//TRIM(VarName)//' not found!')
  END IF

  IF(InfoActive(20)) THEN
    CALL VectorValuesRange(PresSol % Values,SIZE(PresSol % values),PresSol % Name)       
    CALL VectorValuesRange(HeightSol % Values,SIZE(HeightSol % values),HeightSol % Name)       
  END IF

  
  IF( SensMode ) THEN
    VarName = ListGetString( Params,'Sensitivity Variable Name', Found )
    IF(.NOT. Found) THEN
      CALL Fatal(Caller,'Give "Sensitivity Variable Name", e.g. "Displacement Adjoint 3"!')
    END IF
    SensSol => VariableGet( Mesh % Variables, VarName )
    IF(.NOT. ASSOCIATED(SensSol)) THEN
      CALL Fatal(Caller,'"Sensitivity" variable '//TRIM(VarName)//' not found!')
    END IF
    IF(InfoActive(20)) THEN
      CALL VectorValuesRange(SensSol % Values,SIZE(SensSol % values),SensSol % Name)       
    END IF
  END IF
     
  
  grav = 0.0_dp
  IF( ListGetLogical( Params,'Use Gravity',Found ) ) THEN  
    gWork => ListGetConstRealArray( Model % Constants,'Gravity',Found)
    IF(Found) THEN
      grav = ABS(gWork(SIZE(gWork,1),1))
    END IF
  END IF

  i = ListGetInteger( Params,'Fix AC After Iters',Found )
  IF( i > 0 ) THEN
    Var => VariableGet( Mesh % Variables,'coupled iter' )
    IF( NINT(Var % Values(1)) > i ) RETURN
  END IF
    


  
  Var => Solver % Variable  
  IF ( .NOT. AllocationsDone ) THEN
    n = Mesh % MaxElementNodes
    m = SIZE(Var % Values)
    ALLOCATE( FORCE( n ), STIFF(n,n), ElemHeight(n), PrevElemHeight(n), ElemSens(n), ElemPres(n), PrevElemPres(n), &
        Density(n), Thickness(n), EqPres(n), Nodes % X(n), Nodes % Y(n), Nodes % Z(n), STAT=istat ) 
    IF ( istat /= 0 ) CALL Fatal(Caller, 'Memory allocation error.' )
    AllocationsDone = .TRUE.
    LSInit = ListGetLogical( Params,'Levelset Initialize',Found )
  ELSE
    LSInit = .FALSE.
  END IF

  Volume = 0.0_dp
  Area = 0.0_dp
  PresInt = 0.0_dp
      
  IF(LSInit) THEN
    VarName = GetString( Params,'Levelset Variable Name', Found )
    LSVar => VariableGet( Mesh % Variables, VarName, UnfoundFatal = .TRUE.) 

    H3Int = 0.0_dp
    LSMin = 0.0_dp
    DO t=1,Solver % NumberOfActiveElements
      Element => GetActiveElement(t)
      n = GetElementNOFNodes()
      CALL LocalLS( Element, n )
    END DO    
    Area = ParallelReduction(Area)
    H3Int = ParallelReduction(H3Int)
    LSMin = ParallelReduction(LSMin,1)
   
    ! Use the above parameters to set the initial values
    CALL SetACInit()
    CALL Info(caller,'Initialized the artificial compressibility field!',Level=6)
    RETURN   
  END IF
  
  ! Get maximum values for dp and dh since the may affect the formulation for compressibility
  dpmax = 0.0_dp
  dhmax = 0.0_dp
  nlift = 0
  DO t=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()
    CALL LocalMaxima( Element, n, dpmax, dhmax, nlift )
  END DO
  dpmax = ParallelReduction( dpmax, 2 ) 
  dhmax = ParallelReduction( dhmax, 2 ) 

  WRITE(Message,'(A,ES12.3)') 'Maximum pressure overhead:', dpmax
  CALL Info(Caller,Message,Level=10)
  WRITE(Message,'(A,ES12.3)') 'Maximum displacement:', dhmax
  CALL Info(Caller,Message,Level=10)
  CALL Info(Caller,'Number of lifted elements: '//I2S(nlift),Level=10)

  ! We will solve this as a PDE only when there is diffusion, even zero.
  ! Otherwise the matrix equation will be fully diagonal. 
  Diff = ListGetCReal( Params,'AC Diffusion',SolvePDE)
  
  ac_min = ListGetCReal( Params,'AC Minimum Value',Found )

  dp0 = ListGetCReal( Params,'AC Pressure Epsilon',Found )
  IF(.NOT. Found) THEN
    val = ListGetCReal( Params,'AC Pressure Epsilon Relative',Found )
    IF(.NOT. Found) val = 0.01_dp
    dp0 = MAX(EPSILON(dp0),val*dpmax)
  END IF
  
!------------------------------------------------------------------------------
! Initialize the system and do the assembly
!------------------------------------------------------------------------------
  CALL DefaultInitialize()
  
  DO t=1,Solver % NumberOfActiveElements
    Element => GetActiveElement(t)
    n = GetElementNOFNodes()
    NodeIndexes => Element % NodeIndexes

    CALL LocalMatrix(  STIFF, FORCE, Element, n )
    CALL DefaultUpdateEquations( STIFF, FORCE )
  END DO

  CALL DefaultFinishBoundaryAssembly()
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()

  CALL DefaultFinish()

  Volume = ParallelReduction(Volume)
  Area = ParallelReduction(Area)
  PresInt = ParallelReduction(PresInt)
  
  CALL VectorValuesRange(ACinst,SIZE(ACinst),Solver % Variable % Name)       

  CALL ListAddConstReal(Model % Simulation,'res: AC Volume',Volume)
  WRITE(Message,'(A,ES12.5)') 'AC Volume: ',Volume
  CALL Info(Caller, Message)

  CALL ListAddConstReal(Model % Simulation,'res: AC Area',Area)
  WRITE(Message,'(A,ES12.5)') 'AC Area: ',Area
  CALL Info(Caller, Message)

  val = 0.0_dp
  IF(PresInt > EPSILON(val)) val =  Volume / PresInt
  CALL ListAddConstReal(Model % Simulation,'res: AC Average',val)
  WRITE(Message,'(A,ES12.5)') 'AC Average: ',val
  CALL Info(Caller, Message)
  
  DoIt = .TRUE.
  val = ListGetCReal( Params,'AC Update Condition',Found ) 
  IF(Found .AND. val < 0.0_dp ) THEN
    CALL Info(Caller,'Not updating the AC field because of condition!',Level=7) 
    DoIt = .FALSE.
  END IF
    
  val = ListGetCReal( Params,'AC Critical Volume',Found ) 
  IF(Found .AND. Volume < val ) THEN
    CALL Info(Caller,'Not updating the AC field because of small volume!',Level=7) 
    DoIt = .FALSE.
  END IF

  val = ListGetCReal( Params,'AC Critical Area',Found ) 
  IF(Found .AND. Area < val ) THEN
    CALL Info(Caller,'Not updating the AC field because of small area!',Level=7) 
    DoIt = .FALSE.
  END IF

  Var => VariableGet( Mesh % Variables,'AC')
  ACave => Var % Values
  
  IF(DoIt) THEN
    CALL Info(Caller,'Updating AC field for FSI coupling!',Level=7) 

    ! Do some asymmetric relaxation.
    ! Larger values are more robust. Hence one could use larger relaxation
    ! factor for growing and smaller factor for decreasing values. 
    Relax = ListGetCReal( Params,'AC Relaxation Factor',Found )
    IF(.NOT. Found) Relax = 1.0_dp
    RelaxM = ListGetCReal( Params,'AC Relaxation Factor Negative',Found )
    IF(.NOT. Found) RelaxM = Relax
    RelaxP = ListGetCReal( Params,'AC Relaxation Factor Positive',Found )
    IF(.NOT. Found) RelaxP = Relax

    IF( nVisited == 0 ) THEN
      RelaxM = 1.0_dp
      RelaxP = 1.0_dp
    END IF
    
    val = ListGetCReal( Params,'AC Multiplier',Found )
    IF(.NOT. Found ) val = 1.0_dp
    
    IF(RelaxM*RelaxP < 1.0_dp - EPSILON(Relax)) THEN
      ! We cab relax the decreasing and increasing suggested AC values differently.
      ! Idea is that increaing values should maybe react quicker. 
      WHERE ( val * ACinst > ACave ) 
        ACave = ACave + RelaxP * (val*ACinst-ACave)
      ELSE WHERE
        ACave = ACave + RelaxM * (val*ACinst-ACave)
      END WHERE
    ELSE
      ACave = val * ACinst
    END IF

    ! Set minimum threshold
    ACave = MAX(ACave, ac_min)
  END IF

  IF( InfoActive(15)) THEN
    CALL VectorValuesRange(ACave,SIZE(ACave),Var % Name)
  END IF
  nVisited = nVisited + 1

  CALL Info(Caller,'All done for now!',Level=12)
  
  
CONTAINS



!------------------------------------------------------------------------------
  SUBROUTINE LocalMaxima(  Element, n, dpmax, dhmax, nlift )
!------------------------------------------------------------------------------
    TYPE(Element_t), TARGET :: Element
    INTEGER :: n
    REAL(KIND=dp) :: dpmax, dhmax
    INTEGER :: nlift
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: p0, p1, dpres, h0, h1, dh, th
    INTEGER :: t, nt
    TYPE(ValueList_t), POINTER :: Material
    
    CALL GetScalarLocalSolution( elemHeight, tStep=0, UVariable=HeightSol)   
    CALL GetScalarLocalSolution( elemPres, tStep=0, UVariable=PresSol)
      
    IF(DiffMode) THEN
      CALL GetScalarLocalSolution( prevElemHeight, tStep=-1, UVariable=HeightSol)
      CALL GetScalarLocalSolution( prevElemPres, tStep=-1, UVariable=PresSol)
    END IF
    
    Material => GetMaterial( Element )

    EqPres(1:n) = GetReal( Material,'Hydrostatic Pressure',GotEqPres)
    IF(.NOT. Found) THEN
      Density(1:n) = GetReal( Material, 'Ice Density')
      Thickness(1:n) = GetReal( Material, 'Thickness')
    END IF

    nt = 0
    DO t=1,n
      IF(FullMode) THEN
        h0 = 0.0_dp
        IF(GotEqPres) THEN
          p0 = EqPres(t)
        ELSE          
          rho = Density(t)
          th = Thickness(t)
          p0 = rho * grav * th
        END IF
      ELSE
        ! In differential mode we need the data from previous timestep. 
        h0 = prevElemHeight(t)
        p0 = prevElemPres(t)
      END IF
      
      ! The "height" field could have an offset which is not seen in the difference: h1-h0
      h1 = elemHeight(t)       
      dh = (h1-h0)
      
      ! The "pressure" field could have an offset which is not seen in the difference: p1-p0
      p1 = elemPres(t)
      dpres = p1-p0

      dhmax = MAX(dhmax,h1)
      dpmax = MAX(dpmax,dpres)

      IF(h1 * dpres > 1.0e-20 ) nt = nt + 1
    END DO

    ! If the whole element is lifted then add the counter
    IF(nt == n) nlift = nlift + 1

    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMaxima
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE LocalLS( Element, n )
!------------------------------------------------------------------------------
    TYPE(Element_t), TARGET :: Element
    INTEGER :: n
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), Thickness(n), elemLS(n), DetJ, S, LS
    LOGICAL :: Stat
    INTEGER :: t
    REAL(KIND=dp) :: h
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: Material
    
    IP = GaussPoints( Element )

    Nodes % x(1:n) = Model % Nodes % x(Element % NodeIndexes)
    Nodes % y(1:n) = Model % Nodes % y(Element % NodeIndexes)
    Nodes % z(1:n) = 0_dp

    CALL GetScalarLocalSolution( elemLS, UVariable=LSVAr)

    Material => GetMaterial( Element )

    Thickness(1:n) = GetReal( Material, 'Thickness')
    
    DO t=1,IP % n      
      stat = ElementInfo( Element, Nodes, IP % u(t), IP % v(t), IP % w(t), &
          DetJ, Basis )

      S = IP % s(t) * detJ
      
      LS = SUM(Basis(1:n) * elemLS(1:n) )      
      IF( LS < 0.0 ) THEN
        Area = Area + s
        h = SUM(Basis(1:n) * Thickness(1:n) )
        H3Int = H3Int + s / h**3
        LSMin = MIN(LSMin, LS)
      END IF
    END DO
    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalLS
!------------------------------------------------------------------------------


  SUBROUTINE SetACInit()

    REAL(KIND=dp) :: Youngs, Pois, D, Reff, LS, h, Heff
    INTEGER :: i,j
    TYPE(ValueList_t), POINTER :: Material
    
    ! Here we takes exactly the 1st material from the elastic body!!!
    DO i=1,Model % NumberOfMaterials
      Material => Model % Materials(i) % Values
      Youngs = ListGetCReal( Material,'Youngs Modulus',Found)
      Pois = ListGetCReal( Material,'Poisson Ratio',Found)
      IF(Found) EXIT
    END DO

    ! Effective thickness
    Heff = (H3Int/Area)**(-1.0/3)
    
    !  flexural rigidity D with Young's modulus, effective thickness, and Poisson's ratio
    D = Youngs * Heff**3 / (12 * (1-Pois**2))

#if 1
    Reff = SQRT(Area / PI)
#else
    Reff = -LSMin
#endif 
           
    DO i=1,Mesh % NumberOfNodes
      j = Solver % Variable % Perm(i)
      IF(j==0) CYCLE

      ! displacement with unit load
      LS = LSVar % Values( LSVar % Perm(i) )
      IF( LS > 0.0_dp ) THEN
        h = 0.0_dp
      ELSE
        ! This assumes circular plate
        h = ( Reff**4 / 64 * D ) * ( 1 - (LS/LSMin)**2)**2
      END IF
      
      ACave(j) = h
    END DO

  END SUBROUTINE SetACInit
  
  
  
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix(  STIFF, FORCE, Element, n )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: STIFF(:,:), FORCE(:)
    TYPE(Element_t), TARGET :: Element
    INTEGER :: n
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Basis(n), dBasisdx(n,3)
    REAL(KIND=dp) :: DetJ, S, p0, p1, dpres, h0, h1, dh, th, ac
    LOGICAL :: Stat
    INTEGER :: t, p, q, dim, k
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: Material
    
    STIFF = 0.0_dp
    FORCE = 0.0_dp
    
    dim = CoordinateSystemDimension()
    IP = GaussPoints( Element )

    Nodes % x(1:n) = Model % Nodes % x(Element % NodeIndexes)
    Nodes % y(1:n) = Model % Nodes % y(Element % NodeIndexes)
    !IF(dim == 3) Nodes % z(1:n) = Model % Nodes % z(Element % NodeIndexes)
    Nodes % z(1:n) = 0_dp

    CALL GetScalarLocalSolution( elemHeight, tStep=0, UVariable=HeightSol)   
    CALL GetScalarLocalSolution( elemPres, tStep=0, UVariable=PresSol)
      
    IF(DiffMode) THEN
      CALL GetScalarLocalSolution( prevElemHeight, tStep=-1, UVariable=HeightSol)
      CALL GetScalarLocalSolution( prevElemPres, tStep=-1, UVariable=PresSol)
    END IF

    IF( SensMode ) THEN
      CALL GetScalarLocalSolution( elemSens, UVariable=SensSol)   
    END IF
      
    
    Material => GetMaterial( Element )

    IF(FullMode) THEN
      EqPres(1:n) = GetReal( Material,'Hydrostatic Pressure',GotEqPres)
      IF(.NOT. Found) THEN
        Density(1:n) = GetReal( Material, 'Ice Density')
        Thickness(1:n) = GetReal( Material, 'Thickness')
      END IF
    END IF
      
        
    DO t=1,IP % n
      
      stat = ElementInfo( Element, Nodes, IP % u(t), IP % v(t), IP % w(t), &
          DetJ, Basis, dBasisdx )
      S = IP % s(t) * detJ


      h1 = SUM(elemHeight(1:n) * Basis(1:n))      
      p1 = SUM(elemPres(1:n) * Basis(1:n))

      IF(GotEqPres) THEN
        p0 = SUM(Basis(1:n) * EqPres(1:n))
      ELSE          
        rho = SUM(Basis(1:n) * Density(1:n))
        th = SUM(Basis(1:n) * Thickness(1:n))
        p0 = rho * grav * th
      END IF

      
      IF( SensMode ) THEN
        dh = h1
        dpres = p1-p0        
        ac = SUM(elemSens(1:n) * Basis(1:n))      
      ELSE IF(FullMode) THEN
        dh = h1
        dpres = p1-p0        
        ac = dh / MAX(dpres, dp0 )        
      ELSE
        ! In differential mode we need the data from previous timestep. 
        dh = h1 - SUM(prevElemHeight(1:n) * Basis(1:n))
        dpres = p1 - SUM(prevElemPres(1:n) * Basis(1:n))        
        ac = dh / MAX(dpres, dp0 )        
      END IF      
        
!------------------------------------------------------------------------------
!      Finally, the elemental matrix & vector
!------------------------------------------------------------------------------       

      IF( SolvePDE ) THEN
        ! diffusion of AC parameter     
        STIFF(1:n,1:n) = STIFF(1:n,1:n) + &
            S * Diff * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )

        DO p=1,n
          DO q=1,n
            STIFF(p,q) = STIFF(p,q) + s * Basis(p) * Basis(q)
          END DO
        END DO
      ELSE
        DO p=1,n
          STIFF(p,p) = STIFF(p,p) + s * Basis(p)
        END DO
      END IF

      ! Assembly with the actual data.
      FORCE(1:n) = FORCE(1:n) + s * ac

      IF( h1 > EPSILON(h1)) THEN
        Volume = Volume + s * dh
        Area = Area + s
        ! ~pV 
        PresInt = PresInt + s * dpres
      END IF
    END DO

    
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE FilmCompressibility
!------------------------------------------------------------------------------
