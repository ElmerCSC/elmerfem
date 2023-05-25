! ******************************************************************************
! *
! *  Authors: Thomas Zwinger, Rupert Gladstone
! *  Email:   
! *  Web: 
! *
! *  Original Date: 
! *   2021/05/25
! *****************************************************************************
!> Solver Weertman2Coulomb
!> 
!> Converts linear Weertman coefficient (e.g. from inversion) to "Coulomb"
!> sliding law parameters
!> 
!> At time of writing, this document contains notes and equations pertinent to
!> the conversion of coefficients:
!> https://photos.app.goo.gl/HRnwnKtqUJ5NeT7TA
!> 
!> The following "Conversion mode" options are available (see above document
!> for more info):
!>
!> "Threshold"
!> A threshold value of the Weertman sliding coefficient is given.  Either
!> side of this threshold one or other of the "Coulomb" coefficients is held 
!> constant while the other is derived.
!> 
!> "Smooth"
!> A Weertman equation is used to calculate the As "Coulomb" coefficient.  This is
!> then scaled toward zero for regions where the Coulomb limit is approached 
!> (currently using a tanh function based on effective pressure).  The C
!> coefficient is then derived as a function of As.
!> 
!> 
!> For the .sif:
!> 
!> Required variables...
!> 
!> Normal Vector
!> Effective Pressure
!> Flow Solution
!> 
!> Weertman Coefficient (variable name as given in solver params)
!> Coulomb C            (variable name as given in solver params)
!> Coulomb As           (variable name as given in solver params)
!> 
!> Required solver parameters...
!> 
!> Weertman Coefficient Input Variable
!> Coulomb C Output Variable
!> Coulomb As Output Variable
!> Conversion mode
!> 
!> Optional solver parameters... (required if conversion mode = threshold)
!> 
!> Threshold Sliding Coefficient
!> Minimum As
!> Default C
!> Default As
!> 
! Note that the Weertman2Coulumb function that was previously in USF_Sliding
! also contained the option to provide "external pressure" and calculate
! "effective pressure" itself using also the stress vector.  Can be copied
! from repo version f9ca191d4e48eaa39d7a1ab7e if needed.

SUBROUTINE Weertman2CoulombSolver( Model,Solver,dt,TransientSimulation )

  USE DefUtils
  USE SolverUtils
  IMPLICIT NONE

  ! intent in
  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp)  :: dt
  LOGICAL        :: TransientSimulation
  
  ! local variables
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Variable_t), POINTER  :: AsVar, CVar, lwcVar, EPvar, FlowVariable, NormalVar
  INTEGER                    :: nn, ii, DIM
  INTEGER, POINTER           :: CoulombPerm(:), lwcPerm(:), EPperm(:), CPerm(:), &
       FlowPerm(:), NormalPerm(:), AsPerm(:)
  REAL(KIND=dp)              :: BetaSwitch, defaultC, defaultAs, normalvelocity(3), &
       tangentialvelocity(3), tangentialvelocitysquared, normal(3), velo(3), MinAs
  REAL(KIND=dp), POINTER     :: CoulombParam(:), lwcValues(:), AsValues(:), EPvalues(:), &
       FlowValues(:), CValues(:), NormalValues(:)
  LOGICAL                    :: GotIt
  CHARACTER(LEN=MAX_NAME_LEN):: FlowSolverName, Cname, Asname, CoulombVarName, LWCname
  CHARACTER(LEN=MAX_NAME_LEN):: ConversionMode

  DIM = CoordinateSystemDimension()

  SolverParams => GetSolverParams()

  !--------------------------------------------------------------------------------------------
  ! The solver variable is just there to contain a mask determining which Coulomb 
  ! coefficient is being converted (the other is given a default value)
  CoulombVarName = Solver % Variable % Name
  CoulombParam => Solver % Variable % Values
  CoulombPerm => Solver % Variable % Perm

  !--------------------------------------------------------------------------------------------
  ! Variable containing the (optimised) linear weertman coefficient
  LWCname = GetString( SolverParams , 'Weertman Coefficient Input Variable', GotIt )
  IF (.NOT.GotIt) &
    CALL FATAL("Weertman2Coulomb",'Keyword "Weertman Coefficient Input Variable" not found')
  lwcVar =>  VariableGet(Model % mesh % Variables,TRIM(LWCname),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(lwcVar)) &
       CALL FATAL("Weertman2Coulomb","Variable "//TRIM(LWCname)//" not found")
  lwcPerm => lwcVar % Perm
  lwcValues => lwcVar % Values

  !--------------------------------------------------------------------------------------------
  ! Variables to contain converted coefficients
  Cname = GetString( SolverParams , 'Coulomb C Output Variable', GotIt )
  IF (.NOT.GotIt) &
    CALL FATAL("Weertman2Coulomb",'Keyword "Coulomb C Output Variable" not found')
  CVar =>  VariableGet(Model % Variables,TRIM(Cname),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(CVar)) &
       CALL FATAL("Weertman2Coulomb","Variable "//TRIM(Cname)//" not found")
  CPerm => CVar % Perm
  CValues => CVar % Values
  !
  Asname = GetString( SolverParams , 'Coulomb As Output Variable', GotIt )
  IF (.NOT.GotIt) &
    CALL FATAL("Weertman2Coulomb",'Keyword "Coulomb As Output Variable" not found')
  AsVar =>  VariableGet(Model % Variables,TRIM(Asname),UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(AsVar)) &
       CALL FATAL("Weertman2Coulomb","Variable "//TRIM(Asname)//" not found")
  AsPerm => AsVar % Perm
  AsValues => AsVar % Values
  

  !--------------------------------------------------------------------------------------------
  ! Effective pressure and normal variables are required
  EPvar => VariableGet(Model % Mesh % Variables, "Effective Pressure", UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(EPvar)) THEN
    CALL FATAL("Weertman2Coulomb","Variable >Effective Pressure< not found")
  END IF
  EPperm => EPvar % Perm
  EPvalues => EPvar % Values
  !
  NormalVar =>  VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=.TRUE.)
  IF (.NOT.ASSOCIATED(NormalVar)) &
       CALL FATAL("Weertman2Coulomb",'Variable "Normal Vector" not found')
  NormalPerm => NormalVar % Perm
  NormalValues => NormalVar % Values

  !--------------------------------------------------------------------------------------------
  ! Solver params (constants) for use in the conversion
  ConversionMode = GetString(SolverParams, "Conversion mode", GotIt)
  IF (.NOT.GotIt) &
       CALL FATAL("Weertman2Coulomb",'No "Conversion mode" found')
  !
  IF (ConversionMode.EQ."Threshold") THEN
     BetaSwitch = GetConstReal(SolverParams, "Threshold Sliding Coefficient", GotIt)
     IF (.NOT.GotIt) &
          CALL FATAL("Weertman2Coulomb",'No "Threshold Sliding Coefficient" found')
     !
     MinAs = GetConstReal(SolverParams, "Minimum As", GotIt)
     IF (.NOT.GotIt) &
          CALL FATAL("Weertman2Coulomb",'No "Minimum As" found')
     !
     defaultC = GetConstReal(SolverParams,"Default C", GotIt)
     IF (.NOT.GotIt) &
          CALL FATAL("Weertman2Coulomb",'No "Default C" found')
     !
     defaultAs = GetConstReal(SolverParams,"Default As", GotIt)
     IF (.NOT.GotIt) &
          CALL FATAL("Weertman2Coulomb",'No "Default As" found')
     !
  END IF
     
  FlowSolverName = GetString( Model % Solver % Values , 'Flow Solver Name', GotIt )    
  IF (.NOT.Gotit) FlowSolverName = 'Flow Solution'
  FlowVariable => VariableGet( Model % Variables, FlowSolverName, UnFoundFatal=.True.)
  IF (.NOT.ASSOCIATED(FlowVariable)) &
       CALL FATAL("Weertman2Coulomb",'Variable "FlowVariable" not found')
  FlowPerm    => FlowVariable % Perm
  FlowValues  => FlowVariable % Values


  ! Loop over all nodes
  DO nn = 1, Solver % Mesh % Nodes % NumberOfNodes
     
     IF (EPperm(nn).eq.0.0_dp) CYCLE

     velo = 0.0_dp
     DO ii=1, DIM     
        normal(ii) = -NormalValues(DIM*(NormalPerm(nn)-1) + ii)      
        velo(ii) = FlowValues( (DIM+1)*(FlowPerm(nn)-1) + ii )
     END DO
     normalvelocity = SUM(normal(1:DIM)*velo(1:DIM))*normal
     tangentialvelocity = velo - normalvelocity
     tangentialvelocitysquared = SUM(tangentialvelocity(1:DIM)*tangentialvelocity(1:DIM))
     
     SELECT CASE(ConversionMode)

     CASE("Threshold")
        CValues(CPerm(nn))   = defaultC
        AsValues(AsPerm(nn))  = defaultAs

        IF (lwcValues(lwcPerm(nn)) > BetaSwitch) THEN
           CValues(CPerm(nn))   = defaultC
           AsValues(AsPerm(nn)) = 1.0_dp/(lwcValues(lwcPerm(nn))**3.0_dp * tangentialvelocitysquared) &
                - SQRT(tangentialvelocitysquared)/EPvalues(EPperm(nn))**3.0_dp
           CoulombParam(CoulombPerm(nn)) = 1.0
           IF (AsValues(AsPerm(nn)).lt.MinAs) THEN
              AsValues(AsPerm(nn)) = MinAs
              ! Should we leave CValues(CPerm(nn)) at its default or recalculate?
           END IF
        ELSE
           AsValues(AsPerm(nn))  = defaultAs
           CValues(CPerm(nn))    = lwcValues(lwcPerm(nn))*SQRT(tangentialvelocitysquared)/EPvalues(EPperm(nn))
           CoulombParam(CoulombPerm(nn)) = -1.0
        END IF

     CASE("Smooth")

        ! As = u_b^(1-n).beta^-n
        ! (assuming n = 3)
        AsValues(AsPerm(nn)) = &
             tangentialvelocitysquared**(-1) * lwcValues(lwcPerm(nn))**(-3.0_dp)
        CoulombParam(CoulombPerm(nn)) = tanh(2.0_dp*EPvalues(EPperm(nn)))
        AsValues(AsPerm(nn)) = CoulombParam(CoulombPerm(nn)) * AsValues(AsPerm(nn)) 

        ! C = u_b.beta.N^(-1).(1 - beta^3.u_b^2.A_s)^(-1/3)   
        CValues(CPerm(nn)) = SQRT(tangentialvelocitysquared) * lwcValues(lwcPerm(nn))           &
             * EPvalues(EPperm(nn))**(-1.0_dp) * (1.0_dp - lwcValues(lwcPerm(nn))**(3.0_dp)     &
             * SQRT(tangentialvelocitysquared)**2.0_dp * AsValues(AsPerm(nn)) )**(-1.0_dp/3.0_dp) 

     CASE DEFAULT
        CALL FATAL("Weertman2Coulomb",'Conversion mode not recognised')

     END SELECT

  END DO

  NULLIFY(SolverParams)
  NULLIFY(AsVar)
  NULLIFY(CVar)
  NULLIFY(lwcVar)
  NULLIFY(EPvar)
  NULLIFY(FlowVariable)
  NULLIFY(NormalVar)
  NULLIFY(CoulombPerm)
  NULLIFY(lwcPerm)
  NULLIFY(EPperm)
  NULLIFY(CPerm)
  NULLIFY(FlowPerm)
  NULLIFY(NormalPerm)
  NULLIFY(AsPerm)
  NULLIFY(CoulombParam)
  NULLIFY(lwcValues)
  NULLIFY(AsValues)
  NULLIFY(EPvalues)
  NULLIFY(FlowValues)
  NULLIFY(CValues)
  NULLIFY(NormalValues)

END SUBROUTINE Weertman2CoulombSolver

