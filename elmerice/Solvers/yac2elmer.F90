SUBROUTINE YAC2Elmer( Model,Solver,dt,TransientSimulation )
  USE DefUtils, ONLY: GetSolverParams, GetMesh, GetNOFActive, &
    DefaultVariableAdd, GetLogical, GetLogical, GetString, MAX_NAME_LEN, &
    VariableGet, ParEnv, variable_on_elements
  USE GeneralUtils, ONLY: I2S
  USE Types, ONLY: Model_t, Solver_t, Mesh_t, Variable_t, ValueList_t, dp
  USE Messages, ONLY: Message, FATAL, INFO, USE_YAC
  USE elmer_coupling, ONLY: coupling_setup, is_root_rank
  USE elmer_ebfm_coupling, ONLY: elmer_ebfm_interface, t_ice_field, smb_field, &
                                 runoff_field, surface_height_field
  ! USE elmer_icon_coupling, ONLY: elmer_icon_interface, clt_field, pr_field

  IMPLICIT NONE

  TYPE(Model_t),  INTENT(IN) :: Model
  TYPE(Solver_t), INTENT(IN) :: Solver
  REAL(KIND=dp),  INTENT(IN) :: dt
  LOGICAL,        INTENT(IN) :: TransientSimulation


  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t), POINTER :: ThisMesh
  CHARACTER(LEN=MAX_NAME_LEN):: SolverName='YAC2Elmer'
  ! parameters to be read in from this solvers section in the sif
  LOGICAL :: couple_to_ebfm, couple_to_icon         ! define which component is coupled to Elmer

  CHARACTER(LEN=1024) ::  config_file, model_tstep, coupling_timestep, grid_crs
  INTEGER :: I, t, ierr, dt_hours, coupling_hours
  INTEGER, POINTER :: t_icePerm(:), smbPerm(:), runoffPerm(:)
  ! INTEGER, POINTER :: cltPerm(:), prPerm(:)  ! ICON is not supported at the moment
  LOGICAL :: Parallel, FirstTime=.TRUE., UnFoundFatal=.TRUE.
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Variable_t), POINTER :: t_iceVar, smbVar, runoffVar, ZsSol
  ! TYPE(Variable_t), POINTER :: cltVar, prVar  ! ICON is not supported at the moment

  LOGICAL        :: Found

  SolverParams => GetSolverParams()

  ! read config file
  config_file = GetString(SolverParams, 'Config File Name',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >Config File Name< found in yac2elmer solver')
  END IF
  ! check if config file actually exists
  INQUIRE(FILE=TRIM(config_file), EXIST=USE_YAC)
  IF (.NOT. USE_YAC) THEN
    CALL FATAL(SolverName,'Coupling config file not found')
  END IF

  ! check if it is coupled to EBFM
  couple_to_ebfm = GetLogical(SolverParams, 'Couple To EBFM',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >Couple To EBFM< found in yac2elmer solver')
  END IF

  ! check if it is coupled to ICON
  couple_to_icon = GetLogical(SolverParams, 'Couple To ICON',  Found )
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >Couple To ICON< found in yac2elmer solver')
  END IF

  ! read coupling timestep as ISO 8601 string
  coupling_timestep = GetString(SolverParams, 'Coupling Time Step', Found)
  IF (.NOT. Found) THEN
     CALL FATAL(SolverName,'No keyword >Coupling Time Step< found in yac2elmer solver')
  END IF
  ! Consistency check: Coupling Time Step must be ISO 8601 duration string like 'PT3H', 'PT24H', etc.
  IF (.NOT. (LEN_TRIM(coupling_timestep) >= 4 .AND. &
       coupling_timestep(1:2) == 'pt' .AND. &
       coupling_timestep(LEN_TRIM(coupling_timestep):LEN_TRIM(coupling_timestep)) == 'h')) THEN
  CALL FATAL(SolverName, &
    "'Coupling Time Step' must be provided as ISO 8601 duration string of the form " // &
    "'PT3H', 'PT24H', etc. Other formats are currently not supported.")
  END IF
  ! Extract the number between 'PT' and 'H' and convert to integer
  READ(coupling_timestep(3:LEN_TRIM(coupling_timestep)-1), *, IOSTAT=ierr) coupling_hours
  IF (ierr /= 0) THEN
     CALL FATAL(SolverName,"Could not parse number of hours from 'Coupling Time Step'")
  END IF

  ! read grid CRS (coordinate reference system)
  grid_crs = GetString(SolverParams, 'Grid CRS', Found)
  IF (.NOT. Found) THEN
     grid_crs = "EPSG:3413"  ! default value
     CALL INFO(SolverName, 'No keyword >Grid CRS< found, using default: EPSG:3413', Level=3)
  END IF

  IF (.NOT. (couple_to_ebfm .OR. couple_to_icon)) THEN
    CALL FATAL(SolverName,'At least one of >Couple To EBFM< or >Couple To ICON< must be TRUE')
  END IF
  ! TODO: remove this check when ICON coupling is implemented
  IF (couple_to_icon) THEN
    CALL FATAL(SolverName,'>Couple To ICON< is currently not supported. Please set to FALSE')
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! retrieve the timestep in hours
  Mesh => Solver % Mesh

  dt_hours = int(dt * 8760)
  write(model_tstep, *) dt_hours
  CALL INFO(SolverName, &
    'ELMER timestep size in hours:' // I2S(dt_hours), Level=3)
  ! Check if coupling_hours is greater than or equal to dt_hours
  IF (coupling_hours < dt_hours) THEN
     CALL FATAL(SolverName,"'Coupling Time Step' must be greater than or equal to Elmer time step size in hours")
  END IF
  ! Check if coupling_hours is an integer multiple of dt_hours
  IF (MOD(coupling_hours, dt_hours) /= 0) THEN
     CALL FATAL(SolverName,"'Coupling Time Step' must be an integer multiple of Elmer time step size in hours")
  END IF
  IF (FirstTime) THEN


    ! get mesh
    ThisMesh => GetMesh(Solver)

    ! check if this is a parallel run
    IF ((ParEnv % PEs <= 1) .AND. ( .NOT. ThisMesh % SingleMesh )) THEN
      CALL FATAL(SolverName,'Only parallel runs can use this solver')
    ELSE
      CALL INFO(SolverName, &
        'Running with ' // I2S(ParEnv % PEs) // ' partitions', Level=3)
    END IF

    CALL coupling_setup(ThisMesh, TRIM(grid_crs), TRIM(ADJUSTL(I2S(coupling_hours))), couple_to_ebfm, couple_to_icon)

    IF (couple_to_ebfm) THEN
      ! setting up Elmer-side variables for receiving YAC variables
      ! this coul dbe replaced by an automatic picking of the names and the DOFs
      ! from the coupling-deifnitions
      ! nodal variable
      ALLOCATE(t_icePerm(Mesh % NumberOfNodes), runoffPerm(GetNOFActive(Solver)), &
          smbPerm(GetNOFActive(Solver)))
      DO i=1,Mesh % NumberOfNodes
        t_icePerm(i) = i
      END DO

      ! This is for a element(cell) variables
      DO t=1,GetNOFActive(Solver)
        smbPerm(t) = t
        runoffPerm(t) = t
      END DO

      ! nodal variables (everything that is not a flux)
      CALL DefaultVariableAdd('T_ice', dofs=1, Perm = t_icePerm)
      !CALL DefaultVariableAdd('runoff', dofs=1, Perm = runoffPerm)

      ! element wise (cell) variable
      CALL DefaultVariableAdd('smb', dofs=1, VariableType = Variable_on_elements, Perm = smbPerm)
      CALL DefaultVariableAdd('runoff', dofs=1, VariableType = Variable_on_elements, Perm = runoffPerm)

      ! get surface elevation from Elmer for the first time step
      ZsSol => VariableGet( Model % Mesh % Variables, "Zs" ,UnFoundFatal=UnFoundFatal)
      DO i=1, Mesh % NumberOfNodes
        surface_height_field(i,1) = ZsSol % Values(ZsSol % Perm(i))
      END DO
    END IF

    IF (couple_to_icon) THEN
      CALL FATAL(SolverName,'ICON coupling not yet implemented')

      ! ALLOCATE(cltPerm(Mesh % NumberOfNodes), prPerm(GetNOFActive(Solver)))
      ! DO i=1,Mesh % NumberOfNodes
      !   cltPerm(i) = i
      ! END DO

      ! DO t=1,GetNOFActive(Solver)
      !   prPerm(t) = t
      ! END DO

      ! CALL DefaultVariableAdd('tas', dofs=1, Perm = cltPerm)

      ! ! element wise (cell) variable
      ! CALL DefaultVariableAdd('pr_snow', dofs=1, VariableType = Variable_on_elements, Perm = prPerm)
    END IF

    FirstTime = .FALSE.

    CALL INFO(SolverName, "YAC coupling setup done", Level=1)
  END IF
!!!!!!!!!! DO WE HAVE TO INITIALIZE WITH EVERY CALL ? !!!!!!!!!!!!!!

  IF (couple_to_ebfm) THEN
      CALL INFO(SolverName, 'BEFORE ELMER EBFM INTERFACE', Level=3)
      ! couple with EBFM
      CALL elmer_ebfm_interface(is_root_rank)
      CALL INFO(SolverName, 'AFTER ELMER EBFM INTERFACE', Level=3)

      t_iceVar => VariableGet( Mesh % Variables, 'T_ice' )
      smbVar => VariableGet( Mesh % Variables, 'smb' )
      runoffVar => VariableGet( Mesh % Variables, 'runoff' )
      ZsSol => VariableGet( Model % Mesh % Variables, "Zs" ,UnFoundFatal=UnFoundFatal)
      CALL INFO(SolverName, 'AFTER GETTING VARIABLES', Level=3)
      IF ((.NOT.ASSOCIATED(t_iceVar)) .OR. (.NOT.ASSOCIATED(smbVar)) .OR. (.NOT.ASSOCIATED(runoffVar))) THEN
        CALL FATAL(SolverName,'Elmer variables not associated')
      END IF

      CALL INFO(SolverName, 'BEFORE WRITING NODAL VALUES', Level=3)
       !write over values for nodes
      DO i=1, Mesh % NumberOfNodes
        t_iceVar % Values(t_iceVar % Perm(i)) = t_ice_field(i,1)
        surface_height_field(i,1) = ZsSol % Values(ZsSol % Perm(i))
      END DO
      CALL INFO(SolverName, 'BEFORE WRITING ELEMENT VALUES', Level=3)
      ! write over values for elements
      DO t=1, GetNOFActive(Solver)
         smbVar  % Values(smbVar %  Perm(t)) = smb_field(t,1)
         runoffVar  % Values(runoffVar %  Perm(t)) = runoff_field(t,1)
      END DO
      !CALL INFO(SolverName, "Test output start", Level=1)
      !PRINT *, "Size of clt_field", SIZE(clt_field, 1),"First entry:", clt_field(1,1)
      !CALL INFO(SolverName, "Test output start", Level=1)
      !CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )
  END IF

  IF (couple_to_icon) THEN
      CALL FATAL(SolverName,'ICON coupling not yet implemented')
      ! TODO: stub implementation for ICON coupling
      ! CALL elmer_icon_interface(is_root_rank)
      ! cltVar => VariableGet( Mesh % Variables, 'tas' )
      ! prVar => VariableGet( Mesh % Variables, 'pr_snow' )
  END IF

  CALL INFO(SolverName,'Coupling step done', Level=1)

END SUBROUTINE YAC2Elmer

