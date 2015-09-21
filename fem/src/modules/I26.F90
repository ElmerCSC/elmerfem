!/******************************************************************************
! *
! *  Module for defining circuits and dynamic equations
! *
! *  Authors: Juha Ruokolainen, Eelis Takala, Antero Arkkio
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 30.11.2012
! *
! *****************************************************************************/

!> \ingroup Solvers
!> \{


!------------------------------------------------------------------------------
!> Initialization for the primary solver: CurrentSource
!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamics_init( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
  TYPE(ValueList_t), POINTER :: Params

  Params => Solver % Values

  ! When we introduce the variables in this way the variables are created
  ! so that they exist when the proper simulation cycle starts.
  ! This also keeps the command file cleaner.
  CALL ListAddString( Params,'Exported Variable 1',&
      '-global Rotor Angle')
  CALL ListAddString( Params,'Exported Variable 2',&
      '-global Rotor Velo')
  Solver % Values => Params

!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamics_init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamics( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver       !< Linear & nonlinear equation solver options
  TYPE(Model_t) :: Model         !< All model information (mesh, materials, BCs, etc...)
  REAL(KIND=dp) :: dt            !< Timestep size for time dependent simulations
  LOGICAL :: TransientSimulation !< Steady state or transient simulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  LOGICAL :: Found, First=.TRUE.

  TYPE(Solver_t), POINTER :: Asolver => Null(), RotMSolver => Null()

  INTEGER :: p,q,bid,ind,nn,nd,nm,ind0
  TYPE(ValueList_t), POINTER :: BF,Params,BC
  INTEGER, POINTER :: PS(:)
  REAL(KIND=dp)::vemf,vphi,vind,A,b,sigma,SumResistance
  TYPE(Matrix_t), POINTER :: CM=>Null()
  INTEGER :: tstep=0, CompId
  TYPE(Variable_t), POINTER :: LagrangeVar, AngVar, VeloVar
  REAL(KIND=dp), TARGET :: torq,imom=0,ang=0._dp,velo=0._dp,scale,sclA,sclB

  REAL(KIND=dp), ALLOCATABLE, SAVE :: Tcoef(:,:,:), RotM(:,:,:)
  REAL(KIND=DP), POINTER, SAVE :: Cwrk(:,:,:)
  TYPE(Mesh_t), POINTER :: Mesh  
  TYPE(Valuelist_t), POINTER :: Material, BodyParams
  TYPE(Bodyarray_t), POINTER :: CircuitVariableBody

  integer :: slen,i,j,k,m,n,istat,BodyId,ColId,RowId,jj,tot_nofcomp
  CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd,CoilType,dofnumber
  REAL(KIND=dp) :: BodyY

  LOGICAL  :: owner, STAT

  TYPE OldCircuitVariable_t
    LOGICAL :: isIvar, isVvar
    INTEGER :: BodyId, valueId, dofs, pdofs
    TYPE(OldComponent_t), POINTER :: Component
    REAL(KIND=dp), ALLOCATABLE :: A(:), B(:), Source(:)
    INTEGER, ALLOCATABLE :: EqVarIds(:)
  END TYPE

  TYPE(OldCircuitVariable_t), POINTER :: Cvar

  TYPE OldComponent_t
    REAL(KIND=dp) :: ElArea, N_j, coilthickness, nofturns
    INTEGER :: BodyId, polord, ElBoundary, nofcnts
    CHARACTER(LEN=MAX_NAME_LEN) :: CoilType
    TYPE(OldCircuitVariable_t), POINTER :: ivar, vvar
  END TYPE

  TYPE(OldComponent_t), POINTER :: Comp

  TYPE OldCircuit_t
    REAL(KIND=dp), ALLOCATABLE :: A(:,:), B(:,:), Area(:)
    INTEGER, ALLOCATABLE :: Body(:), Perm(:)
    LOGICAL :: UsePerm = .FALSE.
    INTEGER :: n, m, n_comp,CvarDofs
    CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: names(:), source(:), sourcetype(:)
    TYPE(OldComponent_t), POINTER :: Components(:)
    TYPE(OldCircuitVariable_t), POINTER :: CircuitVariables(:)
  END TYPE OldCircuit_t

  REAL(KIND=dp), ALLOCATABLE, SAVE :: ip(:)
  LOGICAL, ALLOCATABLE, SAVE :: Adirichlet(:)
  LOGICAL :: dofsdone

  INTEGER, SAVE :: n_Circuits, OldCircuit_tot_n
  TYPE(OldCircuit_t), ALLOCATABLE, SAVE :: Circuits(:)

  INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
  LOGICAL*1, ALLOCATABLE :: Done(:)
  REAL(KIND=dp), POINTER :: Values(:)

  TYPE(Element_t), POINTER :: e, e_p

  TYPE(Variable_t), POINTER, SAVE :: RotMvar

  INTEGER, PARAMETER :: ind1(9) = [1,1,1,2,2,2,3,3,3]
  INTEGER, PARAMETER :: ind2(9) = [1,2,3,1,2,3,1,2,3]

!------------------------------------------------------------------------------
  
  IF (First) THEN
    First = .FALSE.
    
    Mesh => Model % Mesh
    N = Mesh % MaxElementDOFs

    ALLOCATE( Tcoef(3,3,N), RotM(3,3,N), STAT=istat )
    IF ( istat /= 0 ) THEN
      CALL Fatal( 'WhitneyAVSolver', 'Memory allocation error.' )
    END IF

    NULLIFY( Cwrk )

    ! Look for the solver we attach the circuit equations to:
    ! -------------------------------------------------------
    DO i=1,Model % NumberOfSolvers
      Asolver => Model % Solvers(i)
      IF(ListCheckPresent(Asolver % Values,'Export Lagrange Multiplier'))EXIT
    END DO

    ! Look for the RotM solver:
    ! ----------------------
    DO i=1,Model % NumberOfSolvers
      RotMSolver => Model % Solvers(i)
      IF(ListCheckPresent(RotMSolver % Values,'RotM Solver')) EXIT
    END DO

    ! Initialize circuit matrices:
    ! ----------------------------

    CALL Circuits_Init()

    ! Flag constrained A-field DOFs:
    ! ------------------------------
    ALLOCATE(Adirichlet(aSolver % Matrix % NumberOfRows))
    Adirichlet = .FALSE.
    DO i=1,GetNOFBOundaryElements()
      e => GetBoundaryElement(i)

      BC => GetBC()
      IF(.NOT.ASSOCIATED(BC)) CYCLE

      IF( ListCheckPresent( BC, 'A {e}' ) ) THEN
        e_p => e % BoundaryInfo % Left
        IF(.NOT.ASSOCIATED(e_p)) e_p => e % BoundaryInfo % Right
        IF (.NOT.ASSOCIATED(e_p) ) CYCLE

        e => Find_Face(e_p, e )
        IF(ASSOCIATED(e)) &
          Adirichlet(PS(e % EdgeIndexes+aSolver % Mesh % NumberofNodes)) = .TRUE.
      END IF
    END DO

    imom = GetConstReal( GetSimulation(),'Imom', Found) ! interatial moment of the motor
    IF (.NOT. Found) imom = 0._dp
  END IF

  IF (Tstep /= GetTimestep()) THEN
    Tstep = GetTimestep()

    ! Dynamic equations, export rotation angle:
    ! -----------------------------------------
    torq = GetConstReal( GetSimulation(),'res: Air Gap Torque', Found)

    AngVar => DefaultVariableGet( 'Rotor Angle' )
    ! Variable should alreadt exist as it was introduced in the _init section.
    IF(.NOT. ASSOCIATED( AngVar ) ) THEN
      CALL Fatal('CurrentSource','Variable > Rotor Angle < does not exist!')
    END IF
    ang = AngVar % Values(1)

    VeloVar => DefaultVariableGet( 'Rotor Velo' )
    IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
      CALL Fatal('CurrentSource','Variable > Rotor Velo < does not exist!')
    END IF
    velo = VeloVar % Values(1)

    ! Time integration for the angular momentum equation
    !---------------------------------------------------
    IF(imom /= 0) THEN
      velo = velo + dt * (torq-0) / imom
      ang  = ang  + dt * velo
    END IF

    VeloVar % Values(1) = velo
    AngVar % Values(1) = ang

    ! Circuit variable values from previous timestep:
    ! -----------------------------------------------
    ip = 0._dp
    LagrangeVar => VariableGet( Solver % Mesh % Variables,'LagrangeMultiplier')
    IF(ASSOCIATED(LagrangeVar)) THEN
      IF(SIZE(LagrangeVar % Values)>=OldCircuit_tot_n) ip=LagrangeVar % Values(1:OldCircuit_tot_n)
    END IF

    ! Export circuit & dynamic variables for "SaveScalars":
    ! -----------------------------------------------------

    CALL ListAddConstReal(GetSimulation(),'res: time', GetTime())

    DO p=1,n_Circuits
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i)), ip(Cvar % ValueId))

        IF (Cvar % pdofs /= 0 ) THEN
          DO jj = 1, Cvar % pdofs
            write (dofnumber, "(I2)") jj
            CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))&
                                   //' dof '//TRIM(dofnumber), ip(Cvar % ValueId + jj))
          END DO
        END IF

      END DO
    END DO
    CALL ListAddConstReal(GetSimulation(),'res: Speed', velo/(2._dp*pi)*60)
  END IF


  ! Generate values for the circuit matrix entries:
  ! -----------------------------------------------
!  print *, ParEnv % Mype, "Circuits_Apply"
  CALL Circuits_Apply()
!  print *, ParEnv % Mype, "List_toCRSMatrix"
  IF(  CM % Format == MATRIX_LIST ) CALL List_toCRSMatrix(CM)
!  print *, ParEnv % Mype, "DONE"

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE Circuits_Init()
!------------------------------------------------------------------------------
    ! Read Circuit defintions from MATC:
    ! ----------------------------------
    cmd = "Circuits"
    slen = LEN_TRIM(cmd)
    CALL Matc( cmd, name, slen )
    READ(name, *) n_Circuits

    ALLOCATE( Circuits(n_Circuits) )

    ! Read in number of circuit variables and their names for each circuit. 
    ! ---------------------------------------------------------------------
    OldCircuit_tot_n = 0._dp
    tot_nofcomp=0
    DO p=1,n_Circuits
      ! #variables for circuit "p":
      ! ---------------------------
      cmd = 'C.'//TRIM(i2s(p))//'.variables'
      slen = LEN_TRIM(cmd)
      CALL Matc( cmd, name, slen )

      READ(name, *) Circuits(p) % n

      n = Circuits(p) % n
      ALLOCATE( Circuits(p) % body(n), Circuits(p) % names(n) )
      ALLOCATE( Circuits(p) % source(n), Circuits(p) % sourcetype(n) )
      ALLOCATE( Circuits(p) % CircuitVariables(n), Circuits(p) % Perm(n) )
      Circuits(p) % body = 0
      Circuits(p) % names = ' '

     ! Count bodies and circuit variables in MATC
     ! ------------------------------------------
      Circuits(p) % n_comp = 0
      DO i=1,Circuits(p) % n
        cmd = 'C.'//TRIM(i2s(p))//'.name.'//TRIM(i2s(i))
        slen = LEN_TRIM(cmd)
        CALL Matc( cmd, name, slen )

        Circuits(p) % names(i) = name(1:slen)

        Circuits(p) % CircuitVariables(i) % isIvar = .FALSE.
        Circuits(p) % CircuitVariables(i) % isVvar = .FALSE.

        IF(name(1:7) == 'i_body(' .OR. name(1:7) == 'v_body(') THEN
          DO j=8,slen
            IF(name(j:j)==')') EXIT 
          END DO
          READ(name(8:j-1),*) BodyId

          Circuits(p) % CircuitVariables(i) % BodyId = BodyId

          IF ( .NOT. ANY(Circuits(p) % body == BodyID) ) THEN
            Circuits(p) % n_comp = Circuits(p) % n_comp + 1
            Circuits(p) % body(Circuits(p) % n_comp) = BodyId
          END IF

          SELECT CASE (name(1:7))
          CASE('i_body(')
            Circuits(p) % CircuitVariables(i) % isIvar = .TRUE.
          CASE('v_body(')
            Circuits(p) % CircuitVariables(i) % isVvar = .TRUE.
          CASE DEFAULT
            CALL Fatal('Circuits_Init()', 'Circuit variable should be either i_body or v_body!')
          END SELECT
        ELSE
            Circuits(p) % CircuitVariables(i) % isIvar = .FALSE.
            Circuits(p) % CircuitVariables(i) % isVvar = .FALSE.
            Circuits(p) % CircuitVariables(i) % dofs = 1
            Circuits(p) % CircuitVariables(i) % pdofs = 0
            Circuits(p) % CircuitVariables(i) % BodyId = 0
        END IF
      END DO

      ! create components and count the number of dofs in circuit
      ! ---------------------------------------------------------
      ALLOCATE(Circuits(p) % Components(Circuits(p) % n_comp))
      tot_nofcomp=tot_nofcomp+Circuits(p) % n_comp

      Circuits(p) % CvarDofs = 0
      k = 0
      DO CompId=1,Circuits(p) % n_comp
        Comp => Circuits(p) % Components(CompId)
        Comp % nofcnts = 0
        Comp % BodyId = Circuits(p) % body(CompId)
        Comp % ivar => NULL()
        Comp % vvar => NULL()
        
        DO i=1,Circuits(p) % n
          Cvar => Circuits(p) % CircuitVariables(i)

          IF (Comp % BodyId == Cvar % BodyId) THEN
            IF (Cvar % isIvar) THEN
              Comp % ivar => Cvar
              Cvar % Component => Comp
              CYCLE
            ELSE IF (Cvar % isVvar) THEN
              Comp % vvar => Cvar
              Cvar % Component => Comp
              CYCLE
            END IF
          END IF
        END DO

        IF (.NOT. ASSOCIATED(Comp % ivar) ) THEN
          CALL FATAL('Circuits_Init', 'Current Circuit Variable is not found for Body '//TRIM(i2s(BodyId)))
        ELSE IF (.NOT. ASSOCIATED(Comp % vvar) ) THEN
          CALL FATAL('Circuits_Init', 'Voltage Circuit Variable is not found for Body '//TRIM(i2s(BodyId)))
        END IF

        CircuitVariableBody => Model % Bodies (Comp % BodyId)

        BodyParams => CircuitVariableBody % Values
        IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('Circuits_Init', 'Body Parameters not found!')
        Comp % CoilType = GetString(BodyParams, 'Coil Type', Found)
        IF (.NOT. Found) CALL Fatal ('Circuits_Init', 'Coil Type not found!')
        
        SELECT CASE (Comp % CoilType) 
        CASE ('stranded')
          Comp % nofturns = GetConstReal(BodyParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')

          Comp % ElBoundary = GetInteger(BodyParams, 'Electrode Boundary 1', Found)
          IF (.NOT. Found) THEN 
            Comp % ElArea = GetConstReal(BodyParams, 'Electrode Area', Found)
            IF (.NOT. Found) THEN
              CALL Fatal('Circuits_Init','Electrode Boundary 1 or Electrode Area not found!')
            END IF
          ELSE
            ! Compute Electrode Area Automatically:
            ! -------------------------------------
!              DO t=1,GetNOFBoundaryElements()
!              Element => GetBoundaryElement(t)
          END IF

          Comp % N_j = REAL(Comp % nofturns) / Comp % ElArea
          
          ! Stranded coil has current and voltage 
          ! variables (which both have a dof):
          ! ------------------------------------
          Comp % ivar % dofs = 1
          Comp % vvar % dofs = 1
          Comp % ivar % pdofs = 0
          Comp % vvar % pdofs = 0

        CASE ('massive')
          ! Massive coil has current and voltage 
          ! variables (which both have a dof):
          ! ------------------------------------
          Comp % ivar % dofs = 1
          Comp % vvar % dofs = 1
          Comp % ivar % pdofs = 0
          Comp % vvar % pdofs = 0

        CASE ('foil winding')
          Comp % polord = GetInteger(BodyParams, 'Foil Winding Voltage Polynomial Order', Found)
          IF (.NOT. Found) Comp % polord = 2

          ! Foil winding has current and voltage 
          ! variables. Current has one dof and 
          ! voltage has a polynom for describing the 
          ! global voltage. The polynom has 1+"polynom order"
          ! dofs. Thus voltage variable has 1+1+"polynom order"
          ! dofs (V=V0+V1*alpha+V2*alpha^2+..):
          ! dofs:
          ! V, V0, V1, V2, ...
          ! ------------------------------------
          Comp % ivar % dofs = 1
          Comp % ivar % pdofs = 0
          Comp % vvar % dofs = Comp % polord + 2
          ! polynom dofs:
          ! -------------
          Comp % vvar % pdofs = Comp % polord + 1

          Comp % coilthickness = GetConstReal(BodyParams, 'Coil Thickness', Found)
          IF (.NOT. Found) CALL Fatal('Circuits_Init','Coil Thickness not found!')

          Comp % nofturns = GetConstReal(BodyParams, 'Number of Turns', Found)
          IF (.NOT. Found) CALL Fatal('Circuits_Init','Number of Turns not found!')
 
        END SELECT
        CALL AddVariableToCircuit(Circuits(p), Comp % ivar, k)
        CALL AddVariableToCircuit(Circuits(p), Comp % vvar, k)
      END DO

      ! add variables that are not associated to components (bodies)
      DO i=1,Circuits(p) % n
        Cvar => Circuits(p) % CircuitVariables(i)
        IF (Cvar % isIvar .OR. Cvar % isVvar) CYCLE
        CALL AddVariableToCircuit(Circuits(p), Circuits(p) % CircuitVariables(i), k)
      END DO
      k = k + Circuits(p) % CvarDofs
    END DO

    ! Communicate component numbers and variable dofs to bodies:
    ! ----------------------------------------------------------
    DO p=1,n_Circuits
      DO CompId=1,Circuits(p) % n_comp
   
        Comp => Circuits(p) % Components(CompId)   

        CircuitVariableBody => Model % Bodies (Comp % BodyId)

        BodyParams => CircuitVariableBody % Values
        IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('Circuits_Init', 'Body Parameters not found!')

        CALL listAddInteger(BodyParams, 'Component', CompId)
        CALL listAddInteger(BodyParams, 'Circuit Voltage Variable Id', Comp % vvar % valueId)
        CALL listAddInteger(BodyParams, 'Circuit Voltage Variable dofs', Comp % vvar % dofs)
        CALL listAddInteger(BodyParams, 'Circuit Current Variable Id', Comp % ivar % valueId)
        CALL listAddInteger(BodyParams, 'Circuit Current Variable dofs', Comp % ivar % dofs)
        CALL listAddConstReal(BodyParams, 'Stranded Coil N_j', Comp % N_j)
        Model % Bodies(Comp % BodyId) % Values => BodyParams
      END DO
    END DO

    ! Read in the coefficient matrices for the circuit equations:
    ! Ax' + Bx = source:
    ! ------------------------------------------------------------
    DO p=1,n_Circuits
      n = Circuits(p) % n
      ALLOCATE( Circuits(p) % A(n,n), Circuits(p) % B(n,n) )

      CALL matc_get_array('C.'//TRIM(i2s(p))//'.A'//CHAR(0),Circuits(p) % A,n,n)
      CALL matc_get_array('C.'//TRIM(i2s(p))//'.B'//CHAR(0),Circuits(p) % B,n,n)
      
      Circuits(p) % Perm = -1
      DO i=1,n
        cmd = 'nc:C.'//TRIM(i2s(p))//'.perm('//TRIM(i2s(i-1))//')'
        slen = LEN_TRIM(cmd)
        CALL Matc( cmd, name, slen )
        IF (name(6:10) == 'ERROR') CYCLE
        READ(name(1:5),*) Circuits(p) % Perm(i)
      END DO
      
      IF(ALL(Circuits(p) % Perm /= -1)) THEN 
        Circuits(p) % UsePerm = .TRUE.
        CALL Info( 'CircuitsAndDynamics','Found Permutation vector for circuit '//i2s(p), Level=4 )
      END IF
      
      DO i=1,n
        ! Names of the source functions, these functions should be found
        ! in the "Body Force 1" block of the .sif file.
        ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
        ! ---------------------------------------------------------------------------
        cmd = 'nc:C.'//TRIM(i2s(p))//'.source.'//TRIM(i2s(i))
        slen = LEN_TRIM(cmd)
        CALL Matc( cmd, name, slen )
        Circuits(p) % Source(i) = name(1:slen)

        ! Types of the source functions: voltage or current
        ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
        ! ---------------------------------------------------------------------------
        cmd = 'nc:C.'//TRIM(i2s(p))//'.source.'//TRIM(i2s(i))//'.type'
        slen = LEN_TRIM(cmd)
        CALL Matc( cmd, name, slen )
        Circuits(p) % sourcetype(i) = name(1:slen)

        Cvar => Circuits(p) % CircuitVariables(i)
        RowId = Cvar % ValueId

        ALLOCATE(Cvar % A(Circuits(p) % n), &
                 Cvar % B(Circuits(p) % n), &
                 Cvar % Source(Circuits(p) % n))
        Cvar % A = 0._dp
        Cvar % B = 0._dp
        Cvar % Source = 0._dp

        DO j=1,Circuits(p) % n
          IF(Circuits(p) % A(i,j)/=0._dp.OR.Circuits(p) % B(i,j)/=0._dp) THEN
            IF (Circuits(p) % A(i,j)/=0) Cvar % A(j) = Circuits(p) % A(i,j)
            IF (Circuits(p) % B(i,j)/=0) Cvar % B(j) = Circuits(p) % B(i,j)
          END IF
        END DO

      END DO
    END DO

    ! Create CRS matrix strucures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()

    ! Store values of independent variables (e.g. currents & voltages) from
    ! previous timestep here:
    ! ---------------------------------------------------------------------
    ALLOCATE( ip(OldCircuit_tot_n) ); ip = 0._dp
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_Init
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE CreateCircuitVariable (Variable, Component, isIvar, dofs, valueId)
!------------------------------------------------------------------------------
    TYPE(OldCircuitVariable_t), POINTER :: Variable
    LOGICAL :: isIvar, isVvar
    INTEGER, OPTIONAL :: valueId
    INTEGER, OPTIONAL :: dofs
    TYPE(OldComponent_t), POINTER :: Component
    
    IF (.NOT. ASSOCIATED(Variable) ) THEN
      ALLOCATE(Variable)
    END IF

    Variable % isIvar = isIvar
    Variable % isVvar = .NOT. isIvar
    Variable % BodyId = BodyId
    IF (PRESENT(valueId)) THEN
      Variable % valueId = valueId
    ELSE
      Variable % valueId = 0
    END IF
    IF (PRESENT(dofs))  THEN
      Variable % dofs = dofs
    ELSE
      Variable % dofs = 0
    END IF
    Variable % Component => Component
!------------------------------------------------------------------------------
  END SUBROUTINE CreateCircuitVariable
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE AddVariableToCircuit(Circuit, Variable, k)
!------------------------------------------------------------------------------
    TYPE(OldCircuit_t) :: Circuit
    TYPE(OldCircuitVariable_t) :: Variable
    INTEGER :: k

    IF (Circuit % UsePerm) THEN
      Variable % valueId = Circuit % Perm(OldCircuit_tot_n + 1)
    ELSE
      Variable % valueId = OldCircuit_tot_n + 1
    END IF
    
    OldCircuit_tot_n = OldCircuit_tot_n + Variable % dofs
    
!------------------------------------------------------------------------------
  END SUBROUTINE AddVariableToCircuit
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE CompArea(A,Element,nn,nd)
!------------------------------------------------------------------------------
    INTEGER :: nn, nd
    TYPE(Element_t) :: Element

    REAL(KIND=dp) :: Basis(nd), DetJ,A
    INTEGER :: i,t
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    CALL GetElementNodes( Nodes )

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      stat = ElementInfo(Element,Nodes,IP % U(t),IP % V(t),IP % W(t),detJ,Basis)
      A = A + IP % s(t)*detJ
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CompArea
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Circuits_MatrixInit()
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,p,q,n,temp

    ! Initialialize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows

    CM => AllocateMatrix()
    CM % Format = MATRIX_CRS
    ALLOCATE(CM % RHS(nm + OldCircuit_tot_n)); CM % RHS=0._dp

!    cm % format = matrix_list ! xxxx
!    Asolver %  Matrix % AddMatrix => CM
!return ! xxxx

    CM % NumberOfRows = nm + OldCircuit_tot_n
    n = CM % NumberOfRows
    ALLOCATE(Rows(n+1), Cnts(n)); Rows=0; Cnts=0
    ALLOCATE(Done(nm))

    ! COUNT SIZES:
    ! ============

    ! Basic circuit equations...
    ! ---------------------------
    owner = ParEnv % Mype == ParEnv % PEs-1
    IF ( owner ) THEN
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % n
          Cvar => Circuits(p) % CircuitVariables(i)
          RowId = Cvar % ValueId
          DO j=1,Circuits(p) % n
            IF(Cvar % A(j)/=0._dp.OR.Cvar % B(j)/=0._dp) &
               Cnts(RowId + nm) = Cnts(RowId + nm) + 1
          END DO
        END DO
      END DO
    END IF

    ! ... + the terms including reference to @a/@t:
    ! ---------------------------------------------
    DO p=1,n_Circuits
      DO CompId=1,Circuits(p) % n_comp
        Done = .FALSE.
        Comp => Circuits(p) % Components(CompId)
        Cvar => Comp % vvar
        RowId = Cvar % ValueId
        ColId = Cvar % ValueId

        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          Cnts(RowId + nm) = Cnts(RowId + nm) + 1
          IF (owner) Cnts(RowId + nm) = Cnts(RowId + nm) + 1
        CASE('massive')
          Cnts(RowId + nm) = Cnts(RowId + nm) + 1
          IF (owner) Cnts(RowId + nm) = Cnts(RowId + nm) + 1
        CASE('foil winding')
          IF (owner) THEN
            ! V = V0 + V1*alpha + V2*alpha^2 + ...
            Cnts(RowId + nm) = Cnts(RowId + nm) + Cvar % dofs
          END IF
          ! Circuit eqns for the pdofs:
          ! I(Vj) - I = 0
          ! ------------------------------------
          DO j=1, Cvar % pdofs
            Cnts(j + RowId + nm) = Cnts(j + RowId + nm) + Cvar % dofs
          END DO
        END SELECT

!        temp = SUM(Cnts)
!print *, "Active elements", ParEnv % Mype, ":", GetNOFActive()
        DO q=GetNOFActive(),1,-1
          Element => GetActiveElement(q)
          IF (Element % BodyId == Comp % BodyId) THEN
            nn = GetElementNOFNodes(Element)
            nd = GetElementNOFDOFs(Element,ASolver)
            SELECT CASE (Comp % CoilType)
            CASE('stranded')
              CALL CountAndCreate_Vemf(Element,nn,nd,RowId,Cnts,Jsind=ColId)
            CASE('massive')
              IF (.NOT. HasSupport(Element,nn)) CYCLE
              CALL CountAndCreate_Vemf(Element,nn,nd,RowId,Cnts)
            CASE('foil winding')
              IF (.NOT. HasSupport(Element,nn)) CYCLE 
              DO j = 1, Cvar % pdofs
                dofsdone = ( j==Cvar%pdofs )
                CALL CountAndCreate_Vemf(Element,nn,nd,j+RowId,Cnts,dofsdone=dofsdone)
              END DO
            END SELECT
          END IF
        END DO
!        Comp % nofcnts = SUM(Cnts) - temp
!        print *, ParEnv % Mype, "CompId:", CompId, "Comp % nofcnts", Comp % nofcnts
      END DO
    END DO


    ! ALLOCATE CRS STRUCTURES (if need be):
    ! =====================================

    n = SUM(Cnts)

    IF (n<=0) THEN
      CM % NUmberOfRows = 0
      DEALLOCATE(Rows,Cnts,Done,CM); CM=>Null(); RETURN
    END IF

    ALLOCATE(Cols(n+1), Values(n+1))
    Cols = 0; Values = 0._dp

    ! CREATE ROW POINTERS:
    ! ====================

    CM % NumberOfRows = nm + OldCircuit_tot_n
    Rows(1) = 1
    DO i=2,CM % NumberOfRows+1
      Rows(i) = Rows(i-1) + Cnts(i-1)
    END DO

    Cnts = 0

    ! CREATE COLMUNS:
    ! ===============

    ! Basic circuit equations...
    ! ---------------------------
    owner = ParEnv % Mype == ParEnv % PEs-1
    IF ( owner ) THEN
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % n
          Cvar => Circuits(p) % CircuitVariables(i)
          RowId = Cvar % ValueId
          DO j=1,Circuits(p) % n
            IF(Cvar % A(j)/=0._dp .OR. Cvar % B(j)/=0._dp) THEN
              ColId = Circuits(p) % CircuitVariables(j) % ValueId
              Cols(Rows(RowId + nm) + Cnts(RowId + nm)) = ColId + nm
              Cnts(RowId + nm) = Cnts(RowId + nm) + 1
            END IF
          END DO
        END DO
      END DO
    END IF

    ! ... + the terms including reference to @a/@t:
    ! ---------------------------------------------

    DO p=1,n_Circuits
      DO CompId = 1, Circuits(p) % n_comp
        Done = .FALSE.
        Comp => Circuits(p) % Components(CompId)
        Cvar => Comp % vvar
        RowId = Comp % vvar % ValueId
        ColId = Comp % ivar % ValueId


        SELECT CASE (Comp % CoilType)
        CASE('stranded')
          Cols(Rows(RowId + nm) + Cnts(RowId + nm)) = ColId + nm
          Cnts(RowId + nm) = Cnts(RowId + nm) + 1
          IF (owner) THEN
            Cols(Rows(RowId + nm) + Cnts(RowId + nm)) = RowId + nm
            Cnts(RowId + nm) = Cnts(RowId + nm) + 1
          END IF
        CASE('massive')
          IF (owner) THEN
            Cols(Rows(RowId + nm) + Cnts(RowId + nm)) = ColId + nm
            Cnts(RowId + nm) = Cnts(RowId + nm) + 1
          END IF
          Cols(Rows(RowId + nm) + Cnts(RowId + nm)) = RowId + nm
          Cnts(RowId + nm) = Cnts(RowId + nm) + 1
        CASE('foil winding')
          DO j=0, Cvar % pdofs
            IF (owner) THEN
              ! V = V0 + V1*alpha + V2*alpha^2 + ...
              Cols(Rows(RowId + nm) + Cnts(RowId + nm)) = j + RowId + nm
              Cnts(RowId + nm) = Cnts(RowId + nm) + 1
            END IF
            IF (j/=0) THEN
              ! Circuit eqns for the pdofs:
              ! I(Vi) - I = 0
              ! ------------------------------------
                Cols(Rows(j + RowId + nm) + Cnts(j + RowId + nm)) = ColId + nm
                Cnts(j + RowId + nm) = Cnts(j + RowId + nm) + 1
              DO jj = 1, Cvar % pdofs
                Cols(Rows(j + RowId + nm) + Cnts(j + RowId + nm)) = jj + RowId + nm
                Cnts(j + RowId + nm) = Cnts(j + RowId + nm) + 1
              END DO
            END IF
          END DO
        END SELECT

!        temp = SUM(Cnts)
!print *, "Active elements ", ParEnv % Mype, ":", GetNOFActive()
        DO q=GetNOFActive(),1,-1
          Element => GetActiveElement(q)
          IF (Element % BodyId == Comp % BodyId) THEN
            nn = GetElementNOFNodes(Element)
            nd = GetElementNOFDOFs(Element,ASolver)
            SELECT CASE (Comp % CoilType)
            CASE('stranded')
              CALL CountAndCreate_Vemf(Element,nn,nd,RowId,Cnts,indCol=ColId,Cols=Cols,Jsind=ColId)
            CASE('massive')
              IF (.NOT. HasSupport(Element,nn)) CYCLE 
              CALL CountAndCreate_Vemf(Element,nn,nd,RowId,Cnts,Cols=Cols)
            CASE('foil winding')
              IF (.NOT. HasSupport(Element,nn)) CYCLE 
              DO j = 1, Cvar % pdofs
                dofsdone = ( j==Cvar%pdofs )
                CALL CountAndCreate_Vemf(Element,nn,nd,j+RowId,Cnts,Cols=Cols,dofsdone=dofsdone)
              END DO
            END SELECT
          END IF
        END DO
!        Comp % nofcnts = SUM(Cnts) - temp
!        print *, ParEnv % Mype, "CompId:", CompId, "Coil Type:", Comp % CoilType, &
!                 "Comp % BodyId:", Comp % BodyId, "Comp % nofcnts", Comp % nofcnts

      END DO
    END DO

    IF (n /= SUM(Cnts)) THEN
      print *, "Counted Cnts:", n, "Applied Cnts:", SUM(Cnts)
      CALL Fatal('Circuits_MatrixInit', &
                 'There were different amount of matrix elements than was counted')
    END IF

    DEALLOCATE( Cnts, Done )
    CM % Rows => Rows
    CM % Cols => Cols
    CM % Values => Values
    CALL CRS_SortMatrix(CM)
    
    Asolver %  Matrix % AddMatrix => CM
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_MatrixInit
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreate_Vemf(Element,nn,nd,ind,Cnts,indCol,Cols,dofsdone,Jsind)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ind
    OPTIONAL :: Cols
    INTEGER :: Cols(:), Cnts(:)
    INTEGER :: p,i,iCol,j,Indexes(nd),js
    LOGICAL, OPTIONAL :: dofsdone
    INTEGER, OPTIONAL :: Jsind, indCol

    nd = GetElementDOFs(Indexes,Element,ASolver)
    i = ind + nm
    IF (PRESENT(indCol)) iCol = indCol + nm
    IF (PRESENT(Jsind)) js = Jsind + nm
    DO p=nn,nd
      j = Indexes(p)
      IF(.NOT.Done(j)) THEN
        IF (PRESENT(dofsdone)) THEN
          Done(j) = dofsdone
        ELSE
          Done(j) = .TRUE.
        END IF
        j = PS(j)
        IF(PRESENT(Cols)) THEN
          IF (PRESENT(indCol)) THEN
            Cols(Rows(j)+Cnts(j)) = iCol
          ELSE
            Cols(Rows(j)+Cnts(j)) = i
          END IF
          Cols(Rows(i)+Cnts(i)) = j
        END IF
        Cnts(i) = Cnts(i)+1
        Cnts(j) = Cnts(j)+1
        IF(PRESENT(Jsind)) THEN
          IF(PRESENT(Cols)) Cols(Rows(j)+Cnts(j)) = js
          Cnts(j) = Cnts(j)+1
        END IF
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreate_Vemf
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE Circuits_Apply()
!------------------------------------------------------------------------------

    ! Initialialize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows
    IF(.NOT.ASSOCIATED(CM)) RETURN

    CM % RHS = 0._dp
    IF(ASSOCIATED(CM % Values)) CM % Values = 0._dp

    BF => Model % BodyForces(1) % Values

    ! Basic circuit equations...
    ! ---------------------------
    owner = ParEnv % Mype == ParEnv % PEs-1
    IF ( owner ) THEN
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % n
          Cvar => Circuits(p) % CircuitVariables(i)
          RowId = Cvar % ValueId
          vphi=0._dp
          IF (ASSOCIATED(BF).AND.Circuits(p) % Source(i) /= ' ' ) &
          vphi = GetCReal(BF, Circuits(p) % Source(i), Found)
          IF (Found) THEN 
            Cvar % Source(i) = vphi
          END IF
          
          CM % RHS(RowId + nm) = Cvar % Source(i)
          
          DO j=1,Circuits(p) % n

            ColId = Circuits(p) % CircuitVariables(j) % ValueId

            ! A d/dt(x): (x could be voltage or current):
            !--------------------------------------------
            IF(Cvar % A(j) /= 0._dp) THEN
              CALL AddToMatrixElement(CM, RowId+nm, ColId+nm, Cvar % A(j)/dt)
              CM % RHS(RowId+nm) = CM % RHS(RowId+nm) + Cvar % A(j)*ip(RowId)/dt
            END IF
            ! B x:
            ! ------
            IF(Cvar % B(j) /= 0._dp) THEN
              CALL AddToMatrixElement(CM, RowId+nm, ColId+nm, Cvar % B(j))
            END IF
          END DO
        END DO
      END DO
    END IF

    ! ... + the terms including reference to @a/@t + convert currents to current
    ! densities as source for the vector potential:
    ! --------------------------------------------------------------------------
    DO p=1,n_Circuits
      DO CompID = 1, Circuits(p) % n_comp
        Comp => Circuits(p) % Components(CompId)
        Cvar => Comp % vvar
        RowId = Comp % vvar % ValueId
        ColId = Comp % ivar % ValueId
        
        IF ( owner ) THEN
          SELECT CASE (Comp % CoilType)
          CASE('stranded')
            CALL AddToMatrixElement(CM, RowId+nm, RowId+nm, 1._dp)
          CASE('massive')
            CALL AddToMatrixElement(CM, RowId+nm, ColId+nm, -1._dp)
          CASE('foil winding')
            ! Foil Winding voltage: 
            ! V + ...added next... = 0
            ! ----------------------
            CALL AddToMatrixElement(CM, RowId+nm, RowId+nm, 1._dp)
            DO j = 1, Cvar % pdofs 
              ! Foil Winding voltage: 
              !  ... - Nf/Lalpha * int_0^{Lalpha}(V_0+V_1*alpha+V_2*alpha**2+...) = 0
              !          => ... - Nf * (V_0*Lalpha^0 + V_1/2*Lalpha^1 + V_2/3*Lalpha^2 + ...) = 0
              ! where V_m is the mth dof of the polynomial
              ! --------------------------------------------------------------
              CALL AddToMatrixElement(CM, RowId+nm, j + RowId+nm, &
                  -REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1))

              ! Circuit eqns for the pdofs:
              ! - Nf/Lalpha * I * int_0^1(Vi'(alpha)) + ...added later... = 0
              ! ----------------------------------------------------------
              CALL AddToMatrixElement(CM, j+RowId+nm, ColId+nm, &
                 - REAL(Comp % nofturns) / REAL(j) * Comp % coilthickness**(j-1))
            END DO
          END SELECT
        END IF

        DO q=GetNOFActive(),1,-1
          Element => GetActiveElement(q)
          IF (Element % BodyId == Comp % BodyId) THEN
            BodyParams => GetBodyParams( Element )
            IF (.NOT. ASSOCIATED(BodyParams)) CALL Fatal ('Circuits_apply', 'Body Parameters not found')

            CoilType = GetString(BodyParams, 'Coil Type', Found)
            IF (.NOT. Found) CoilType = ''

            nn = GetElementNOFNodes(Element)
            nd = GetElementNOFDOFs(Element,ASolver)
            CALL GetConductivity()
            SELECT CASE(CoilType)
            CASE ('stranded')
              CALL Add_stranded(Element,Tcoef(1,1,1:nn),Comp,nn,nd)
            CASE ('massive')
              IF (.NOT. HasSupport(Element,nn)) CYCLE
              CALL Add_massive(Element,Tcoef(1,1,1:nn),Comp,nn,nd)
            CASE ('foil winding')
              IF (.NOT. HasSupport(Element,nn)) CYCLE 
              CALL GetElementRotM(Element, RotM, nn)
              CALL Add_foil_winding(Element,Tcoef,RotM,Comp,nn,nd)
            CASE DEFAULT
              CALL Fatal ('Circuits_apply', 'Non existent Coil Type Chosen!')
            END SELECT
          END IF
        END DO
      END DO
    END DO
    Asolver %  Matrix % AddMatrix => CM

!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_Apply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE GetConductivity()
!------------------------------------------------------------------------------

    Tcoef = 0.0d0
    Material => GetMaterial( Element )
    IF (.NOT. ASSOCIATED(Material)) CALL Fatal('Circuits_apply','Material not found.')

    CALL ListGetRealArray( Material, &
           'Electric Conductivity', Cwrk, nn, Element % NodeIndexes, Found )

    IF (.NOT. Found) CALL Fatal('Circuits_apply', 'Electric Conductivity not found.')

    IF (Found) THEN
       IF ( SIZE(Cwrk,1) == 1 ) THEN
          DO i=1,3
             Tcoef( i,i,1:nn ) = Cwrk( 1,1,1:nn )
          END DO
       ELSE IF ( SIZE(Cwrk,2) == 1 ) THEN
          DO i=1,MIN(3,SIZE(Cwrk,1))
             Tcoef(i,i,1:nn) = Cwrk(i,1,1:nn)
          END DO
       ELSE
          DO i=1,MIN(3,SIZE(Cwrk,1))
             DO j=1,MIN(3,SIZE(Cwrk,2))
                Tcoef( i,j,1:nn ) = Cwrk(i,j,1:nn)
             END DO
          END DO
       END IF
    END IF


!------------------------------------------------------------------------------
  END SUBROUTINE GetConductivity
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Add_stranded(Element,Tcoef,Comp,nn,nd)
!------------------------------------------------------------------------------
    INTEGER :: nn, nd
    REAL(KIND=dp) :: Tcoef(nn)
    TYPE(Element_t) :: Element
    TYPE(OldComponent_t) :: Comp

    REAL(KIND=dp) :: Basis(nn), DetJ,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nn,3), WBasis(nd,3), RotWBasis(nd,3), &
                     wBase(nn), w(3), localC
    INTEGER :: p,q,i,j,t,Indexes(nd),RowId,ColId
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(pPOT,'A',Element,ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,'A',Element,ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF

    CALL GetLocalSolution(Wbase,'W')

    RowId = Comp % vvar % ValueId + nm
    ColId = Comp % ivar % ValueId + nm

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      w = -MATMUL(WBase(1:nn), dBasisdx(1:nn,:))
      w = w/SQRT(SUM(w**2._dp))
      localC = SUM(Tcoef(1:nn) * Basis(1:nn))

      CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)

      ! I * R, where 
      ! R = (1/sigma * js,js):
      ! ----------------------

      CALL AddToMatrixElement(CM, RowId, ColId, &
            Comp % N_j * IP % s(t)*detJ*SUM(w*w)/localC)
      
      DO p=nn+1,nd
        q=p-nn

        IF (Comp % N_j/=0._dp) THEN
          ! ( a,w )
            CALL AddToMatrixElement(CM, RowId, PS(Indexes(p)), &
                   tscl * Comp % N_j * IP % s(t)*detJ*SUM(WBasis(q,:)*w)/localC)
            CM % RHS(RowId) = CM % RHS(RowId) + &
                  Comp % N_j * IP % s(t)*detJ*pPOT(p)*SUM(WBasis(q,:)*w)/localC 

          IF (.NOT. Adirichlet(PS(indexes(p)))) THEN
            ! source: 
            ! (J, rot a'), where
            ! J = w*I, thus I*(w, rot a'):
            ! ----------------------------
            CALL AddToMatrixElement(CM,PS(Indexes(p)), ColId, &
               -Comp % N_j*IP % s(t)*detJ*SUM(WBasis(q,:)*w))
          END IF 
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_stranded
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_massive(Element,Tcoef,Comp,nn,nd)
!------------------------------------------------------------------------------
    INTEGER :: nn, nd
    REAL(KIND=dp) :: Tcoef(nn)
    TYPE(Element_t) :: Element
    TYPE(OldComponent_t) :: Comp

    REAL(KIND=dp) :: Basis(nn), DetJ,POT(nd),pPOT(nd),ppPOT(nd),tscl, localC,gradv(3)
    REAL(KIND=dp) :: dBasisdx(nn,3), WBasis(nd,3), RotWBasis(nd,3),wBase(nn),w(3)
    INTEGER :: p,q,i,j,t,Indexes(nd),pp,RowId
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(pPOT,'A',Element,ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,'A',Element,ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF

    CALL GetLocalSolution(Wbase,'W')

    RowId = Comp % vvar % ValueId

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      localC = SUM(Tcoef(1:nn) * Basis(1:nn))

      CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)

      gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))

      ! computing the source term Vi(sigma grad v0, grad si):
      ! ------------------------------------------------

      CALL AddToMatrixElement(CM, RowId+nm, RowId+nm, &
              IP % s(t)*detJ*localC*SUM(gradv*gradv))
      
      DO j=1,nd-nn
        q = j+nn
        
        CALL AddToMatrixElement(CM, RowId+nm, PS(Indexes(q)), &
               tscl * IP % s(t)*detJ*localC*SUM(Wbasis(j,:)*gradv)/dt)

         CM % RHS(RowId+nm) = CM % RHS(RowId+nm) &
                + IP % s(t)*detJ*localC*pPOT(q)*SUM(WBasis(j,:)*gradv)/dt

        IF (.NOT. Adirichlet(PS(indexes(q)))) THEN
          ! computing the source term Vi(sigma grad v0, a'):
          ! ------------------------------------------------
          CALL AddToMatrixElement(CM, PS(indexes(q)), RowId+nm, &
                IP % s(t)*detJ*localC*SUM(gradv*Wbasis(j,:)))
        END IF

      END DO
    END DO

!------------------------------------------------------------------------------
   END SUBROUTINE Add_massive
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE Add_foil_winding(Element,Tcoef,RotM,Comp,nn,nd)
!------------------------------------------------------------------------------
    INTEGER :: nn, nd
    REAL(KIND=dp) :: Tcoef(:,:,:), C(3,3)
    TYPE(Element_t) :: Element
    TYPE(OldComponent_t) :: Comp

    REAL(KIND=dp) :: Basis(nn), DetJ,POT(nd),pPOT(nd),ppPOT(nd),tscl, &
                     localAlpha, localV, localVtest, gradv(3)
    REAL(KIND=dp) :: dBasisdx(nn,3), WBasis(nd,3), RotWBasis(nd,3), & 
                     wBase(nn),w(3),alpha(nn),RotM(3,3,nn),RotMLoc(3,3)
    INTEGER :: p,q,i,j,t,Indexes(nd),pp,RowId,vpolord_tot, &
               vpolord, vpolordtest, dofId, dofIdtest
    LOGICAL :: stat
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)
    CALL GetLocalSolution(pPOT,'A',Element,ASolver,tstep=-1)

    IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
      tscl=1.0_dp
    ELSE
      tscl=1.5_dp
      CALL GetLocalSolution(ppPOT,'A',Element,ASolver,tstep=-2)
      pPot = 2*pPOT - 0.5_dp*ppPOT
    END IF

    CALL GetLocalSolution(Wbase,'W')
    CALL GetLocalSolution(alpha,'Alpha')

    RowId = Comp % vvar % ValueId
    vpolord_tot = Comp % vvar % pdofs - 1

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)

      gradv = MATMUL( WBase(1:nn), dBasisdx(1:nn,:))
      localAlpha = SUM(alpha(1:nn) * Basis(1:nn))
      
      ! alpha is normalized to be in [0,1] thus, 
      ! it needs to be multiplied by the thickness of the coil 
      ! to get the real alpha:
      ! ------------------------------------------------------
      localAlpha = localAlpha * Comp % coilthickness

      ! Compute the conductivity tensor
      ! -------------------------------
      DO i=1,3
        DO j=1,3
          C(i,j) = SUM( Tcoef(i,j,1:nn) * Basis(1:nn) )
          RotMLoc(i,j) = SUM( RotM(i,j,1:nn) * Basis(1:nn) )
        END DO
      END DO

      ! Transform the conductivity tensor:
      ! ----------------------------------
      C = MATMUL(MATMUL(RotMLoc, C),TRANSPOSE(RotMLoc))

      DO vpolordtest=0,vpolord_tot ! V'(alpha)

        localVtest = localAlpha**vpolordtest
        dofIdtest = vpolordtest + 1 + RowId

        DO vpolord = 0, vpolord_tot ! V(alpha)

          localV = localAlpha**vpolord
          dofId = vpolord + 1 + RowId

          ! Computing the stiff term (sigma V(alpha) grad v0, V'(alpha) grad si):
          ! ---------------------------------------------------------------------
          CALL AddToMatrixElement(CM, dofIdtest+nm, dofId+nm, &
                IP % s(t)*detJ*localV*localVtest*SUM(MATMUL(C,gradv)*gradv))
        END DO

        DO j=1,nd-nn
          q = j+nn
          ! computing the mass term (sigma d/dt(a), V'(alpha) grad si):
          ! ---------------------------------------------------------
          CALL AddToMatrixElement(CM, dofIdtest+nm, PS(Indexes(q)), &
               tscl * IP % s(t)*detJ*localVtest*SUM(MATMUL(C,Wbasis(j,:))*gradv)/dt)

          CM % RHS(dofIdtest+nm) = CM % RHS(dofIdtest+nm) &
               + IP % s(t)*detJ*localVtest*pPOT(q)*SUM(MATMUL(C,WBasis(j,:))*gradv)/dt
        END DO

      END DO

      DO vpolord = 0, vpolord_tot ! V(alpha)

        localV = localAlpha**vpolord
        dofId = vpolord + 1 + RowId

        DO j=1,nd-nn
          q = j+nn
          IF (.NOT. Adirichlet(PS(indexes(q)))) THEN
            ! computing the source term (sigma V(alpha) grad v0, a'):
            ! -------------------------------------------------------
            CALL AddToMatrixElement(CM, PS(indexes(q)), dofId+nm, &
                    IP % s(t)*detJ*localV*SUM(MATMUL(C,gradv)*Wbasis(j,:)))
          END IF
        END DO
      END DO

    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_foil_winding
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE GetElementRotM(Element,RotM,n)
!------------------------------------------------------------------------------
   TYPE(Element_t) :: Element
   INTEGER :: k, l, m, j, n
   REAL(KIND=dp) :: RotM(3,3,n)
   TYPE(Variable_t), POINTER, SAVE :: RotMvar
   LOGICAL, SAVE :: visited = .FALSE.
 

   IF(.NOT. visited) THEN
     visited = .TRUE.
     RotMvar => VariableGet( Mesh % Variables, 'RotM E')
     IF(.NOT. ASSOCIATED(RotMVar)) THEN
       CALL Fatal('GetElementRotM','RotM E variable not found')
     END IF
   END IF

   RotM = 0._dp
   DO j = 1, n
     DO k=1,RotMvar % DOFs
       RotM(ind1(k),ind2(k),j) = RotMvar % Values( &
             RotMvar % DOFs*(RotMvar % Perm(Element % DGIndexes(j))-1)+k)
     END DO
   END DO

!------------------------------------------------------------------------------
 END SUBROUTINE GetElementRotM
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   FUNCTION HasSupport(Element, nn) RESULT(support)
!------------------------------------------------------------------------------
    INTEGER :: nn
    TYPE(Element_t) :: Element
    TYPE(OldComponent_t) :: Comp
    LOGICAL :: support
    REAL(KIND=dp) :: wBase(nn)
    
    CALL GetLocalSolution(Wbase,'W')
    
    support = .FALSE. 
    IF ( ANY(Wbase .ne. 0d0) ) support = .TRUE.
!------------------------------------------------------------------------------
   END FUNCTION HasSupport
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamics
!------------------------------------------------------------------------------
