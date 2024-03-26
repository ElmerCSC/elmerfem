!/******************************************************************************
! *
! *  Module for defining circuits and dynamic equations.
! *  This is the original version of circuit simulator still used by some.
! *
! *  Authors: Juha Ruokolainen
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
!> Initialization for the primary solver
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
  TYPE(Solver_t), POINTER :: Asolver => Null()
  INTEGER :: p,q,bid,ind,nn,nd,nm,ind0
  TYPE(ValueList_t), POINTER :: BF,Params
  INTEGER, POINTER :: PS(:)
  REAL(KIND=dp)::vemf,vphi,vind,A,b,sigma, ang0=0._dp, velo0=0._dp
  TYPE(Matrix_t), POINTER :: CM=>Null()
  INTEGER :: tstep=0
  TYPE(Variable_t), POINTER, SAVE :: LagrangeVar, AngVar, VeloVar
  REAL(KIND=dp), TARGET :: torq,imom=0,ang=0._dp,velo=0._dp,scale,sclA,sclB

  INTEGER :: slen,i,j,k,m,n
  CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd

  LOGICAL  :: owner

  TYPE(ValueList_t), POINTER :: BC
  TYPE(Element_t), POINTER :: e,e_p
  LOGICAL, ALLOCATABLE, SAVE :: Adirichlet(:)

  TYPE Circuit_tt
    REAL(KIND=dp), ALLOCATABLE :: A(:,:), B(:,:), Area(:)
    INTEGER, ALLOCATABLE :: Body(:)
    INTEGER :: n, m
    CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: names(:), source(:)
  END TYPE Circuit_tt

  REAL(KIND=dp), ALLOCATABLE, SAVE :: ip(:)

  INTEGER, SAVE :: n_Circuits, Circuit_tot_m
  TYPE(Circuit_tt), ALLOCATABLE, SAVE :: Circuits(:)

  INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
  LOGICAL*1, ALLOCATABLE :: Done(:)
  REAL(KIND=dp), POINTER :: Values(:)
  
  CHARACTER(LEN=MAX_NAME_LEN) :: sname
!------------------------------------------------------------------------------


  IF (First) THEN
    First = .FALSE.

    ! Look for the solver we attach the circuit equations to:
    DO i=1,Model % NumberOfSolvers      
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      j = INDEX( sname,'MagnetoDynamics2DHarmonic')
      IF( j > 0 ) CYCLE
      k = INDEX( sname,'MagnetoDynamics2D')
      IF( k > 0 ) THEN
        ASolver => Model % Solvers(i) 
        EXIT
      END IF
    END DO
    IF(.NOT. ASSOCIATED( ASolver ) ) THEN    
      DO i=1,Model % NumberOfSolvers
        Asolver => Model % Solvers(i)
        IF(ListCheckPresent(Asolver % Values,'Export Lagrange Multiplier'))EXIT
      END DO
    END IF    
    CALL Info('Circuits2D','Associated circuit with solver index: '//I2S(i),Level=10)
    
    AngVar => DefaultVariableGet( 'Rotor Angle' )
    ! Variable should already exist as it was introduced in the _init section.
    IF(.NOT. ASSOCIATED( AngVar ) ) THEN
      CALL Fatal('CurrentSource','Variable > Rotor Angle < does not exist!')
    END IF
    ang = AngVar % Values(1)

    VeloVar => DefaultVariableGet( 'Rotor Velo' )
    IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
      CALL Fatal('CurrentSource','Variable > Rotor Velo < does not exist!')
    END IF
    velo = VeloVar % Values(1)

    ! Initialize circuit matrices:
    ! ----------------------------
    CALL Circuits_Init()
    imom = GetConstReal( GetSimulation(),'Imom') ! interatial moment of the motor
  END IF

  IF (Tstep /= GetTimestep()) THEN
    Tstep = GetTimestep()

    AngVar => DefaultVariableGet( 'Rotor Angle' )
    ! Variable should already exist as it was introduced in the _init section.
    IF(.NOT. ASSOCIATED( AngVar ) ) THEN
      CALL Fatal('CurrentSource','Variable > Rotor Angle < does not exist!')
    END IF
    ang0 = AngVar % Values(1)
    AngVar % Values(1) = ang

    VeloVar => DefaultVariableGet( 'Rotor Velo' )
    IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
      CALL Fatal('CurrentSource','Variable > Rotor Velo < does not exist!')
    END IF
    velo0 = VeloVar % Values(1)
    VeloVar % Values(1) = velo

    ! Circuit variable values from previous timestep:
    ! -----------------------------------------------
    ip = 0._dp
    
    
    sname = LagrangeMultiplierName(ASolver)
    LagrangeVar => VariableGet( ASolver % Mesh % Variables, sname, ThisOnly = .TRUE.)    
    IF(ASSOCIATED(LagrangeVar)) THEN
      IF(SIZE(LagrangeVar % Values)>=Circuit_tot_m) ip=LagrangeVar % Values(1:Circuit_tot_m)
    END IF
   
    ! Export circuit & dynamic variables for "SaveScalars":
    ! -----------------------------------------------------
    CALL ListAddConstReal(GetSimulation(),'res: time', GetTime())
    j = 0
    DO p=1,n_Circuits
      DO i=1,Circuits(p) % m
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i)), ip(i+j))
      END DO
      j = j + Circuits(p) % m
    END DO
    CALL ListAddConstReal(GetSimulation(),'res: Angle(rad)', ang)
    CALL ListAddConstReal(GetSimulation(),'res: Speed(rpm)', velo/(2._dp*pi)*60)
  END IF

  ! Time integration for the angular momentum equation
  !---------------------------------------------------
  torq = GetConstReal( GetSimulation(),'res: Air Gap Torque', Found)

  IF(imom /= 0) THEN
    velo = velo0 + dt * (torq-0) / imom
    ang  = ang0  + dt * velo
  END IF

  VeloVar % Values(1) = velo
  AngVar % Values(1) = ang

  ! Generate values for the circuit matrix entries:
  ! -----------------------------------------------
  CALL Circuits_Apply()
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE Circuits_Init()
!------------------------------------------------------------------------------
    real(kind=dp) :: xxx

    ! Read Circuit definitions from MATC:
    ! ----------------------------------
    slen =  Matc("Circuits", name)
    READ(name(1:slen), *) n_Circuits

    ALLOCATE( Circuits(n_Circuits) )

    ! Read in number of circuit variables and their names for each circuit. Variables may
    ! include references to "u_emf(body)"'s which are special kind of variables attached to
    ! volume integral of time derivative of the vector potential over a "body". These variables
    ! don't have independent equations, instead the corresponding matrix coefficients on the
    ! rows of "u_emf" give the conversion factors from current -> current density over the 
    ! "body" (modulo dividing by volume, which is done here).
    ! ---------------------------------------------------------------------------------------
    Circuit_tot_m = 0._dp
    DO p=1,n_Circuits
      ! #variables for circuit "p":
      ! ---------------------------
      slen = Matc('C.'//i2s(p)//'.variables', name)

      READ(name(1:slen), *) Circuits(p) % n

      n = Circuits(p) % n
      ALLOCATE( Circuits(p) % body(n), Circuits(p) % names(n) )
      ALLOCATE( Circuits(p) % source(n) )
      Circuits(p) % body = 0
      Circuits(p) % names = ' '

      ! names of the variables:
     ! -----------------------
      DO i=1,Circuits(p) % n
        slen = Matc('C.'//i2s(p)//'.name.'//i2s(i),name)

        Circuits(p) % names(i) = name(1:slen)

        IF(name(1:6) == 'u_emf(') THEN
          DO j=7,slen
            IF(name(j:j)==')') EXIT 
          END DO
          READ(name(7:j-1),*) Circuits(p) % body(i)
        END IF
      END DO

      ! #non-u_emf variables in a circuit...
      ! ------------------------------------
      Circuits(p) % m = COUNT(Circuits(p) % body==0)

      ! #non-u_emf variables overall:
      ! ----------------------------
      Circuit_tot_m = Circuit_tot_m + Circuits(p) % m
    END DO

    ! Read in the coefficient matrices for the circuit equations:
    ! Ax' + Bx = source:
    ! ------------------------------------------------------------
    DO p=1,n_Circuits
      n = Circuits(p) % n
      m = Circuits(p) % m
      ALLOCATE( Circuits(p) % A(n,n), Circuits(p) % B(n,n) )

      CALL matc_get_array('C.'//i2s(p)//'.A'//CHAR(0),Circuits(p) % A,n,n)
      CALL matc_get_array('C.'//i2s(p)//'.B'//CHAR(0),Circuits(p) % B,n,n)

      DO i=1,n
        ! Names of the source functions, these functions should be found
        ! in the "Body Force 1" block of the .sif file.
        ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
        ! ---------------------------------------------------------------------------
        slen = Matc('nc:C.'//i2s(p)//'.source.'//i2s(i),name)
        Circuits(p) % Source(i) = name(1:slen)
      END DO
    END DO

    ! Compute the volume each current is attached to:
    ! -----------------------------------------------
    DO p=1,n_Circuits ! loop over circuits
!      print *, '*****  Circuit index: p =', p

      ALLOCATE(Circuits(p) % Area(Circuits(p) % m))
      Circuits(p) % Area = 0

      DO i=1,Circuits(p) % m   ! loop over currents
!        print *, '*****    Current index: i =', i

        DO j=Circuits(p) % m+1,Circuits(p) % n ! loop over bodies
!          print *, '*****      Body index: j =', j
          IF(Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE
!         print *,'*****        Body number: cirquits(i) % body(j) =', circuits(p) % body(j)

          DO q=1,GetNOFActive()
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes()
              nd = GetElementNOFDOFs()
              CALL CompArea(Circuits(p) % Area(i),Element,nn,nd)
            END IF
         END DO ! q
        END DO ! j

        Circuits(p) % Area(i) = ParallelReduction(Circuits(p) % Area(i))
!       print *, '*****    Area: cirquits(p) % area(i) =', Circuits(p) % Area(i)

     END DO ! i
  END DO ! p

    ! Create CRS matrix structures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()

    ! Store values of independent variables (e.g. currents & voltages) from
    ! previous timestep here:
    ! ---------------------------------------------------------------------
    ALLOCATE( ip(Circuit_tot_m) ); ip = 0._dp
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_Init
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
    INTEGER :: i,j,k,p,q,n

    ! Initialilize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows

    CM => AllocateMatrix()
    CM % Format = MATRIX_CRS
    ALLOCATE(CM % RHS(nm + Circuit_tot_m)); CM % RHS=0._dp

    CM % NumberOfRows = nm + Circuit_tot_m
    n = CM % NumberOfRows
    ALLOCATE(Rows(n+1), Cnts(n)); Rows=0; Cnts=0
    ALLOCATE(Done(nm))

    ! COUNT SIZES:
    ! ============

    ! Basic circuit equations...
    ! ---------------------------
    owner = ( ParEnv % Mype == ParEnv % PEs-1 )
    IF ( owner ) THEN
      k = 0
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % m
          DO j=1,Circuits(p) % m
            IF(Circuits(p) % A(i,j)/=0._dp.OR.Circuits(p) % B(i,j)/=0._dp) THEN
              Cnts(i+k+nm) = Cnts(i+k+nm)+1
            END IF
          END DO
        END DO
        k = k + Circuits(p) % m
      END DO
    END IF

    ! ... + the terms including reference to @a/@t:
    ! ---------------------------------------------
    k = 0
    DO p=1,n_Circuits
      DO i=1,Circuits(p) % m
        Done = .FALSE.
        DO j=Circuits(p) % m+1,Circuits(p) % n
          IF (Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE

          DO q=GetNOFActive(),1,-1
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes(Element)
              nd = GetElementNOFDOFs(Element,ASolver)
              CALL CountAndCreate_Vemf(Element,nn,nd,i+k,Cnts)
            END IF
          END DO
        END DO
      END DO
      k = k + Circuits(p) % m
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

    CM % NumberOfRows = nm + Circuit_tot_m
    Rows(1) = 1
    DO i=2,CM % NumberOfRows+1
      Rows(i) = Rows(i-1) + Cnts(i-1)
    END DO
    Cnts = 0

    ! CREATE COLUMNS:
    ! ===============

    ! Basic circuit equations...
    ! ---------------------------
    owner = ( ParEnv % Mype == ParEnv % PEs-1 )
    IF ( owner ) THEN
      k = 0
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % m

          DO j=1,Circuits(p) % m
            IF(Circuits(p) % A(i,j)/=0._dp.OR.Circuits(p) % B(i,j)/=0._dp) THEN
              Cols(Rows(i+k+nm)+Cnts(i+k+nm)) = j+k+nm
              Cnts(i+k+nm) = Cnts(i+k+nm) + 1
            END IF
          END DO
        END DO
        k = k + Circuits(p) % m
      END DO
    END IF

    ! ... + the terms including reference to @a/@t:
    ! ---------------------------------------------

    k = 0
    DO p=1,n_Circuits
      DO i=1,Circuits(p) % m
        Done = .FALSE.
        DO j=Circuits(p) % m+1,Circuits(p) % n
          IF (Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE

          DO q=GetNOFActive(),1,-1
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes(Element)
              nd = GetElementNOFDOFs(Element,ASolver)
              CALL CountAndCreate_Vemf(Element,nn,nd,i+k,Cnts,Cols)
            END IF
          END DO
        END DO
      END DO
      k = k + Circuits(p) % m
    END DO

    DEALLOCATE( Cnts, Done )
    CM % Rows => Rows
    CM % Cols => Cols
    CM % Values => Values
    CALL CRS_SortMatrix(CM)
    
    Asolver %  Matrix % AddMatrix => CM
  END SUBROUTINE Circuits_MatrixInit
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE CountAndCreate_Vemf(Element,nn,nd,ind,Cnts,Cols)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ind, n_start
    OPTIONAL :: Cols
    INTEGER :: Cols(:), Cnts(:)
    TYPE(Mesh_t), POINTER :: Mesh

    INTEGER :: p,i,j,Indexes(nd)

    nd = GetElementDOFs(Indexes,Element,ASolver)

    Mesh => GetMesh()
    n_start = 0
    IF (ALL(Indexes(1:nd)>Mesh % NumberOfNodes)) n_start=0

    i = ind + nm
    DO p=n_start+1,nd
      j = PS(Indexes(p))
      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        IF(PRESENT(Cols)) THEN
          Cols(Rows(j)+Cnts(j)) = i
          Cols(Rows(i)+Cnts(i)) = j
        END IF
        Cnts(i) = Cnts(i)+1
        Cnts(j) = Cnts(j)+1
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreate_Vemf
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Circuits_Apply()
!------------------------------------------------------------------------------
    ! Initialilize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows
    IF(.NOT.ASSOCIATED(CM)) RETURN

    CM % RHS = 0._dp
    CM % Values = 0._dp

    ! Basic circuit equations...
    ! ---------------------------
    owner = ( ParEnv % Mype == ParEnv % PEs-1 )
    IF ( owner ) THEN
      k = 0
      BF => Model % BodyForces(1) % Values
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % m
          vphi = 0._dp
          IF (ASSOCIATED(BF).AND.Circuits(p) % Source(i) /= ' ' ) &
            vphi = GetCReal(BF, Circuits(p) % Source(i), Found)

          CM % RHS(i+k+nm) = dt*vphi
          DO j=1,Circuits(p) % m

            IF(TransientSimulation.AND.Circuits(p) % A(i,j)/=0._dp) THEN
              CALL AddToMatrixElement(CM, i+k+nm, j+k+nm, Circuits(p) % A(i,j))
              CM % RHS(i+k+nm) = CM % RHS(i+k+nm) + Circuits(p) % A(i,j)*ip(i+k)
            END IF

            IF(Circuits(p) % B(i,j) /= 0._dp) THEN
              CALL AddToMatrixElement(CM, i+k+nm, j+k+nm, dt*Circuits(p) % B(i,j))
            END IF
          END DO
        END DO
        k = k + Circuits(p) % m
      END DO
    END IF

    ! ... + the terms including reference to @a/@t + convert currents to current
    ! densities as source for the vector potential:
    ! --------------------------------------------------------------------------
    k = 0
    DO p=1,n_Circuits ! loop over circuits
!     print *, '+++++  Circuit index: p =', p

      DO i=1,Circuits(p) % m ! loop over currents
!        print *, '+++++    Current index: i =', i

        A = Circuits(p) % Area(i)
        DO j=Circuits(p) % m+1,Circuits(p) % n ! loop over bodies
!          print *, '+++++      Body index: j =', j

          IF (Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE
!         print *,'+++++        Body number: cirquits(i) % body(j) =', circuits(p) % body(j)

          sclA = Circuits(p) % B(i,j) / A
          sclB = Circuits(p) % B(j,i) / A

          DO q=GetNOFActive(),1,-1
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes(Element)
              nd = GetElementNOFDOFs(Element,ASolver)
              CALL Add_vemf(Element,nn,nd,i+k,sclA,sclB)
            END IF
          END DO
        END DO
      END DO
      k = k + Circuits(p) % m
    END DO
    Asolver %  Matrix % AddMatrix => CM
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_Apply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Add_vemf(Element,nn,nd,ind,N_u,N_j)
!------------------------------------------------------------------------------
    INTEGER :: nn, nd, ind
    REAL(KIND=dp) :: N_u, N_j
    TYPE(Element_t) :: Element

    REAL(KIND=dp) :: Basis(nn), DetJ,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nn,3), WBasis(nd,3), RotWBasis(nd,3),wBase(nd),w(3)
    LOGICAL :: stat, Wfound, PiolaT, Found
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    INTEGER :: p,q,i,j,t,Indexes(nd),n_start=0

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    Mesh => GetMesh()

    IF(TransientSimulation) THEN
      CALL GetLocalSolution(pPOT,UElement=Element,USolver=ASolver,tstep=-1)

      IF(Solver % Order<2.OR.GetTimeStep()<=2) THEN 
        tscl=1.0_dp
      ELSE
        tscl=1.5_dp
        CALL GetLocalSolution(ppPOT,UElement=Element,USolver=ASolver,tstep=-2)
        pPot = 2*pPOT - 0.5_dp*ppPOT
      END IF
    END IF

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------

      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
                IP % W(t), detJ, Basis,dBasisdx )

      DO q=1,nd

        IF (TransientSimulation) THEN
          CALL AddToMatrixElement(CM, ind+nm, PS(Indexes(q)), &
                 tscl * N_u* IP % s(t)*detJ*Basis(q) )

          CM % RHS(ind+nm) = CM % RHS(ind+nm) + &
                N_u * IP % s(t)*detJ*pPOT(q)*Basis(q)
        END IF

        IF(N_j/=0) THEN
          CALL AddToMatrixElement(CM,PS(Indexes(q)), ind+nm, &
             -N_j*IP % s(t)*detJ*Basis(q) )
        END IF
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_vemf
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE CircuitsAndDynamics



!------------------------------------------------------------------------------
!> Initialization for the primary solver
!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamicsHarmonic_init( Model,Solver,dt,TransientSimulation )
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
END SUBROUTINE CircuitsAndDynamicsHarmonic_init
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
SUBROUTINE CircuitsAndDynamicsHarmonic( Model,Solver,dt,TransientSimulation )
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

  TYPE(Solver_t), POINTER :: Asolver => Null()

  INTEGER :: p,q,bid,ind,nn,nd,nm,ind0
  TYPE(ValueList_t), POINTER :: BF,Params
  INTEGER, POINTER :: PS(:)
  REAL(KIND=dp)::vemf,vphi_re,vphi_im,vind,A,b,sigma
  TYPE(Matrix_t), POINTER :: CM=>Null()
  INTEGER :: tstep=0
  TYPE(Variable_t), POINTER :: LagrangeVar, AngVar, VeloVar
  REAL(KIND=dp), TARGET :: torq,imom=0,ang=0._dp,velo=0._dp,scale,sclA,sclB, omega=1

  integer :: slen,i,j,k,m,n
  CHARACTER(LEN=MAX_NAME_LEN) :: name,cmd

  LOGICAL  :: owner

  TYPE Circuit_tt
    REAL(KIND=dp), ALLOCATABLE :: A(:,:), B(:,:), Area(:)
    REAL(KIND=dp) :: Omega
    LOGICAL :: FoundOmega
    INTEGER, ALLOCATABLE :: Body(:)
    INTEGER :: n, m
    CHARACTER(LEN=MAX_NAME_LEN), ALLOCATABLE :: names(:), source(:)
  END TYPE Circuit_tt

  REAL(KIND=dp), ALLOCATABLE, SAVE :: ip(:)

  LOGICAL :: PiolaVersion = .FALSE.

  INTEGER, SAVE :: n_Circuits, Circuit_tot_m
  TYPE(Circuit_tt), ALLOCATABLE, SAVE :: Circuits(:)

  INTEGER, POINTER :: Rows(:), Cols(:), Cnts(:)
  LOGICAL*1, ALLOCATABLE :: Done(:)
  REAL(KIND=dp), POINTER :: Values(:)

  CHARACTER(LEN=MAX_NAME_LEN) :: sname
  !------------------------------------------------------------------------------

  CALL Info('CircuitsAndDynamics','-------------------------------------------',Level=8 )
  CALL Info('CircuitsAndDynamics','Assembling electric circuit equations',Level=5 )
  
  omega = ListGetConstReal(Model % Simulation, 'Supply Angular Frequency', found)
  IF(.NOT.Found) omega = GetAngularFrequency(Found=Found)

  IF(.NOT.found) THEN
     PRINT *,'**************************************************************'
     PRINT *,'* Supply Angular Frequency was not found in Simulation block *'
     PRINT *,'**************************************************************'
     CALL Fatal('CurrentSource','Variable > Supply Angular Frequency < does not exist!')
  END if

  IF (First) THEN
    First = .FALSE.
    NULLIFY(ASolver)
    
    ! Look for the solver we attach the circuit equations to:
    ! -------------------------------------------------------
    DO i=1,Model % NumberOfSolvers      
      sname = GetString(Model % Solvers(i) % Values, 'Procedure', Found)
      j = INDEX( sname,'MagnetoDynamics2DHarmonic')
      IF( j > 0 ) THEN
        ASolver => Model % Solvers(i) 
        EXIT
      END IF
    END DO
    IF(.NOT. ASSOCIATED( ASolver ) ) THEN    
      DO i=1,Model % NumberOfSolvers
        Asolver => Model % Solvers(i)
        IF(ListCheckPresent(Asolver % Values,'Export Lagrange Multiplier'))EXIT
      END DO
    END IF

    CALL Info('HarmonicCircuits2D','Associated circuit with solver index: '//I2S(i),Level=10)
    
    PiolaVersion = GetLogical(Asolver % Values, 'Use Piola Transform',Found)

    ! Initialize circuit matrices:
    ! ----------------------------
    CALL Circuits_Init()
  END IF

  ip = 0._dp
  sname = LagrangeMultiplierName( ASolver )
  LagrangeVar => VariableGet( ASolver % Mesh % Variables,sname,ThisOnly=.TRUE.)
  IF(ASSOCIATED(LagrangeVar)) THEN
    IF(SIZE(LagrangeVar % Values)>=2*Circuit_tot_m) ip=LagrangeVar % Values(1:2*Circuit_tot_m)
  END IF

  ! Export circuit & dynamic variables for "SaveScalars":
  ! -----------------------------------------------------
  owner = ( Parenv % Mype == ParEnv % PEs-1 )
  CALL ListAddConstReal(GetSimulation(),'res: time', GetTime())
  j = 0
  DO p=1,n_Circuits
    DO i=1,Circuits(p) % m
      IF(owner) THEN
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))//'_re', ip(2*(i+j-1)+1))
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))//'_im', ip(2*(i+j-1)+2))
      ELSE
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))//'_re', 0._dp )
        CALL ListAddConstReal( GetSimulation(), 'res: '//TRIM(Circuits(p) % names(i))//'_im', 0._dp )
      END IF
    END DO
    j = j + Circuits(p) % m
  END DO

  ! Generate values for the circuit matrix entries:
  ! -----------------------------------------------
  CALL Circuits_Apply()
  IF (ASSOCIATED(CM)) THEN
    IF( CM % Format == MATRIX_LIST) CALL List_toCRSMatrix(CM)
  END IF

CONTAINS

!------------------------------------------------------------------------------
  SUBROUTINE Circuits_Init()
!------------------------------------------------------------------------------
    ! Read Circuit definitions from MATC:
    ! ----------------------------------
    slen = Matc("Circuits", name)
    READ(name(1:slen), *) n_Circuits

    ALLOCATE( Circuits(n_Circuits) )

    ! Read in number of circuit variables and their names for each circuit. Variables may
    ! include references to "u_emf(body)"'s which are special kind of variables attached to
    ! volume integral of time derivative of the vector potential over a "body". These variables
    ! don't have independent equations, instead the corresponding matrix coefficients on the
    ! rows of "u_emf" give the conversion factors from current -> current density over the 
    ! "body" (modulo dividing by volume, which is done here).
    ! ---------------------------------------------------------------------------------------
    Circuit_tot_m = 0._dp
    DO p=1,n_Circuits
      ! #variables for circuit "p":
      ! ---------------------------
      slen = Matc('C.'//i2s(p)//'.variables',name)
      READ(name(1:slen), *) Circuits(p) % n

      n = Circuits(p) % n
      ALLOCATE( Circuits(p) % body(n), Circuits(p) % names(n) )
      ALLOCATE( Circuits(p) % source(n) )
      Circuits(p) % body = 0
      Circuits(p) % names = ' '

      ! names of the variables:
     ! -----------------------
      DO i=1,Circuits(p) % n
        slen = Matc('C.'//i2s(p)//'.name.'//i2s(i),name)
        Circuits(p) % names(i) = name(1:slen)

        IF(name(1:6) == 'u_emf(') THEN
          DO j=7,slen
            IF(name(j:j)==')') EXIT 
          END DO
          READ(name(7:j-1),*) Circuits(p) % body(i)
        END IF
      END DO

      ! #non-u_emf variables in a circuit...
      ! ------------------------------------
      Circuits(p) % m = COUNT(Circuits(p) % body==0)

      ! #non-u_emf variables overall:
      ! ----------------------------
      Circuit_tot_m = Circuit_tot_m + Circuits(p) % m
    END DO

    ! Read in the coefficient matrices for the circuit equations:
    ! Ax' + Bx = source:
    ! ------------------------------------------------------------
    DO p=1,n_Circuits
      n = Circuits(p) % n
      m = Circuits(p) % m
      ALLOCATE( Circuits(p) % A(n,n), Circuits(p) % B(n,n) )

      CALL matc_get_array('C.'//i2s(p)//'.A'//CHAR(0),Circuits(p) % A,n,n)
      CALL matc_get_array('C.'//i2s(p)//'.B'//CHAR(0),Circuits(p) % B,n,n)

      DO i=1,n
        ! Names of the source functions, these functions should be found
        ! in the "Body Force 1" block of the .sif file.
        ! (nc: is for 'no check' e.g. don't abort if the MATC variable is not found!)
        ! ---------------------------------------------------------------------------
        slen = Matc('nc:C.'//i2s(p)//'.source.'//i2s(i),name)
        Circuits(p) % Source(i) = name(1:slen)
      END DO
    END DO

    ! Get circuit angular frequency, either from the 
    ! circuit defs, or from the sif-file defs:
    ! -----------------------------------------------
    DO p=1,n_Circuits
      slen = Matc('nc:C.'//i2s(p)//'.omega',name)
      Circuits(p) % FoundOmega = slen >= 1
      IF(Circuits(p) % FoundOmega) &
          READ( name(1:slen), *) Circuits(p) % Omega

      Found = .FALSE.
      DO i=1,Circuits(p) % m
        DO j=Circuits(p) % m+1,Circuits(p) % n
          IF(Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE
          DO q=1,GetNOFActive()
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              IF(Circuits(p) % FoundOmega) THEN
                CALL ListAddConstReal( GetEquation(), &
                      'Angular Frequency', Circuits(p) % Omega )
              ELSE
                Circuits(p) % Omega = GetAngularFrequency( &
                  Found = Circuits(p) % FoundOmega, UElement = Element )
              END IF
              Found = .TRUE.
              EXIT
            END IF
          END DO
          IF( Found) EXIT
        END DO
        IF( Found) EXIT
      END DO
    END  DO

    ! Compute the volume each current is attached to:
    ! -----------------------------------------------
    DO p=1,n_Circuits
      ALLOCATE(Circuits(p) % Area(Circuits(p) % m))
      Circuits(p) % Area = 0

      DO i=1,Circuits(p) % m
        DO j=Circuits(p) % m+1,Circuits(p) % n
          IF(Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE
          DO q=1,GetNOFActive()
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes()
              nd = GetElementNOFDOFs()
              CALL CompArea(Circuits(p) % Area(i),Element,nn,nd)
            END IF
          END DO
        END DO
        Circuits(p) % Area(i) = ParallelReduction(Circuits(p) % Area(i))
      END DO
    END DO

    ! Create CRS matrix structures for the circuit equations:
    ! ------------------------------------------------------
    CALL Circuits_MatrixInit()

    ! Store values of independent variables (e.g. currents & voltages) from
    ! previous timestep here:
    ! ---------------------------------------------------------------------
    ALLOCATE( ip(2*Circuit_tot_m) ); ip = 0._dp
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_Init
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
    INTEGER :: i,j,k,p,q,n

    ! Initialilize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows

    CM => AllocateMatrix()
    CM % Format = MATRIX_LIST
    ALLOCATE(CM % RHS(nm + 2*Circuit_tot_m)); CM % RHS=0._dp

return

    CM % NumberOfRows = nm + Circuit_tot_m
    n = CM % NumberOfRows
    ALLOCATE(Rows(n+1), Cnts(n)); Rows=0; Cnts=0
    ALLOCATE(Done(nm))

    ! COUNT SIZES:
    ! ============

    ! Basic circuit equations...
    ! ---------------------------
    owner = ( ParEnv % Mype == ParEnv % PEs-1 )
    IF ( owner ) THEN
      k = 0
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % m
          DO j=1,Circuits(p) % m
            IF(Circuits(p) % A(i,j)/=0._dp.OR.Circuits(p) % B(i,j)/=0._dp) THEN
              Cnts(i+k+nm) = Cnts(i+k+nm)+1
            END IF
          END DO
        END DO
        k = k + Circuits(p) % m
      END DO
    END IF

    ! ... + the terms including reference to @a/@t:
    ! ---------------------------------------------
    k = 0
    DO p=1,n_Circuits
      Done = .FALSE.
      DO i=1,Circuits(p) % m
        DO j=Circuits(p) % m+1,Circuits(p) % n
          IF (Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE

          DO q=GetNOFActive(),1,-1
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes(Element)
              nd = GetElementNOFDOFs(Element,ASolver)
              CALL CountAndCreate_Vemf(Element,nn,nd,i+k,Cnts)
            END IF
          END DO
        END DO
      END DO
      k = k + Circuits(p) % m
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

    CM % NumberOfRows = nm + Circuit_tot_m
    Rows(1) = 1
    DO i=2,CM % NumberOfRows+1
      Rows(i) = Rows(i-1) + Cnts(i-1)
    END DO
    Cnts = 0

    ! CREATE COLMUNS:
    ! ===============

    ! Basic circuit equations...
    ! ---------------------------
    owner = ( ParEnv % Mype == ParEnv % PEs-1 )
    IF ( owner ) THEN
      k = 0
      DO p = 1,n_Circuits
        DO i=1,Circuits(p) % m

          DO j=1,Circuits(p) % m
            IF(Circuits(p) % A(i,j)/=0._dp.OR.Circuits(p) % B(i,j)/=0._dp) THEN
              Cols(Rows(i+k+nm)+Cnts(i+k+nm)) = j+k+nm
              Cnts(i+k+nm) = Cnts(i+k+nm) + 1
            END IF
          END DO
        END DO
        k = k + Circuits(p) % m
      END DO
    END IF

    ! ... + the terms including reference to @a/@t:
    ! ---------------------------------------------

    k = 0
    DO p=1,n_Circuits
      Done = .FALSE.
      DO i=1,Circuits(p) % m
        DO j=Circuits(p) % m+1,Circuits(p) % n
          IF (Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE

          DO q=GetNOFActive(),1,-1
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes(Element)
              nd = GetElementNOFDOFs(Element,ASolver)
              CALL CountAndCreate_Vemf(Element,nn,nd,i+k,Cnts,Cols)
            END IF
          END DO
        END DO
      END DO
      k = k + Circuits(p) % m
    END DO

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
   SUBROUTINE CountAndCreate_Vemf(Element,nn,nd,ind,Cnts,Cols)
!------------------------------------------------------------------------------
    TYPE(Element_t) :: Element
    INTEGER :: nn, nd, ind
    OPTIONAL :: Cols
    INTEGER :: Cols(:), Cnts(:)

    INTEGER :: p,i,j,Indexes(nd)

    nd = GetElementDOFs(Indexes,Element,ASolver)
    i = ind + nm
    DO p=1,nd
      j = Indexes(p)
      IF(.NOT.Done(j)) THEN
        Done(j) = .TRUE.
        j = PS(j)
        IF(PRESENT(Cols)) THEN
          Cols(Rows(j)+Cnts(j)) = i
          Cols(Rows(i)+Cnts(i)) = j
        END IF
        Cnts(i) = Cnts(i)+1
        Cnts(j) = Cnts(j)+1
      END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE CountAndCreate_Vemf
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE Circuits_Apply()
!------------------------------------------------------------------------------
    integer :: row, col

    ! Initialilize Circuit matrix:
    ! -----------------------------
    PS => Asolver % Variable % Perm
    nm =  Asolver % Matrix % NumberOfRows
    IF(.NOT.ASSOCIATED(CM)) RETURN

    IF(CM % Format==MATRIX_CRS) THEN
      IF(.NOT.ASSOCIATED(CM % Values)) RETURN
      CM % Values = 0._dp
    END IF
    CM % RHS = 0._dp

    ! Basic circuit equations...
    ! ---------------------------
    owner = ( ParEnv % Mype == ParEnv % PEs-1 )
    IF ( owner ) THEN
      k = 0
      BF => Model % BodyForces(1) % Values
      DO p = 1,n_Circuits
        IF(Circuits(p) % FoundOmega) omega = Circuits(p) % omega
        WRITE(Message, * ) 'Circuit(' // I2S(p) // ') Angular Frequency', Omega
        CALL Info( 'ApplyCircuits', Message, Level = 8 )
        DO i=1,Circuits(p) % m
          vphi_re = 0._dp
          vphi_im = 0._dp
          IF (ASSOCIATED(BF).AND.Circuits(p) % Source(i) /= ' ' ) THEN
            vphi_re = GetCReal(BF, TRIM(Circuits(p) % Source(i))//' re', Found)
            vphi_im = GetCReal(BF, TRIM(Circuits(p) % Source(i))//' im', Found)
          END IF

          row = nm + 2*(i + k - 1)
          CM % RHS(row+1) = vphi_re
          CM % RHS(row+2) = vphi_im
          DO j=1,Circuits(p) % m
            col = nm + 2*(j + k - 1)

            IF(Circuits(p) % A(i,j) /= 0._dp) THEN
              CALL AddToMatrixElement(CM, row+1, col+1, 0._dp )
              CALL AddToMatrixElement(CM, row+2, col+2, 0._dp )
              CALL AddToMatrixElement(CM, row+1, col+2, -omega*Circuits(p) % A(i,j) )
              CALL AddToMatrixElement(CM, row+2, col+1,  omega*Circuits(p) % A(i,j) )
            END IF

            IF(Circuits(p) % B(i,j) /= 0._dp) THEN
              CALL AddToMatrixElement(CM, row+1, col+1, Circuits(p) % B(i,j))
              CALL AddToMatrixElement(CM, row+2, col+2, Circuits(p) % B(i,j))
              CALL AddToMatrixElement(CM, row+1, col+2, 0._dp )
              CALL AddToMatrixElement(CM, row+2, col+1, 0._dp )
            END IF
          END DO
        END DO
        k = k + Circuits(p) % m
      END DO
    END IF

    ! ... + the terms including reference to @a/@t + convert currents to current
    ! densities as source for the vector potential:
    ! --------------------------------------------------------------------------
    k = 0
    DO p=1,n_Circuits
      IF(Circuits(p) % FoundOmega) omega = Circuits(p) % omega
      DO i=1,Circuits(p) % m
        A = Circuits(p) % Area(i)
        DO j=Circuits(p) % m+1,Circuits(p) % n
          IF (Circuits(p) % B(i,j)==0._dp.AND.Circuits(p) % B(j,i)==0._dp) CYCLE

          sclA = Circuits(p) % B(i,j) / A * Omega
          sclB = Circuits(p) % B(j,i) / A

          DO q=GetNOFActive(),1,-1
            Element => GetActiveElement(q)
            IF (Element % BodyId == Circuits(p) % body(j)) THEN
              nn = GetElementNOFNodes(Element)
              nd = GetElementNOFDOFs(Element,ASolver)
              CALL Add_vemf(Element,nn,nd,i+k,sclA,sclB)
            END IF
          END DO
        END DO
      END DO
      k = k + Circuits(p) % m
    END DO
    Asolver %  Matrix % AddMatrix => CM
!------------------------------------------------------------------------------
  END SUBROUTINE Circuits_Apply
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE Add_vemf(Element,nn,nd,ind,N_u,N_j)
!------------------------------------------------------------------------------
    INTEGER :: nn, nd, ind
    REAL(KIND=dp) :: N_u, N_j
    TYPE(Element_t) :: Element

    REAL(KIND=dp) :: Basis(nn), DetJ,POT(nd),pPOT(nd),ppPOT(nd),tscl
    REAL(KIND=dp) :: dBasisdx(nn,3), WBasis(nd,3), RotWBasis(nd,3),wBase(nd),w(3)
    INTEGER :: p,q,i,j,t,Indexes(nd), row, col
    LOGICAL :: stat, Wfound
    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP

    CALL GetElementNodes(Nodes)
    nd = GetElementDOFs(Indexes,Element,ASolver)

    CALL GetLocalSolution(Wbase,'W')
    Wfound = .TRUE.
    IF(ALL(Wbase(1:nn)==0._dp)) Wfound=.FALSE.
    W = [0._dp, 0._dp, 1._dp]

    ! Numerical integration:
    ! ----------------------
    IP = GaussPoints(Element)
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis )

      row = 2*(ind-1)+nm
      DO p=1,nd

        col = 2*(PS(Indexes(p))-1)
        CALL AddToMatrixElement(CM, row+1, col+1, 0._dp )
        CALL AddToMatrixElement(CM, row+2, col+2, 0._dp )
        CALL AddToMatrixElement(CM, row+1, col+2, -N_u* IP % s(t)*detJ*Basis(p))
        CALL AddToMatrixElement(CM, row+2, col+1,  N_u* IP % s(t)*detJ*Basis(p))

        CALL AddToMatrixElement(CM, col+1, row+1, -N_j*IP % s(t)*detJ*Basis(p))
        CALL AddToMatrixElement(CM, col+2, row+2, -N_j*IP % s(t)*detJ*Basis(p))
        CALL AddToMatrixElement(CM, col+1, row+2, 0._dp )
        CALL AddToMatrixElement(CM, col+2, row+1, 0._dp )
      END DO
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE Add_vemf
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE CircuitsAndDynamicsHarmonic
!------------------------------------------------------------------------------
