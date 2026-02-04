SUBROUTINE YAC2Elmer( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE SolverUtils
  USE elmer_coupling, ONLY: coupling_setup
  USE elmer_ebfm_coupling, ONLY: elmer_ebfm_interface, t_ice_field, smb_field, &
                                 runoff_field, surface_height_field
  
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp)  :: dt
  LOGICAL        :: TransientSimulation

  
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t), POINTER :: ThisMesh
  CHARACTER(LEN=MAX_NAME_LEN):: SolverName='YAC2Elmer'
  ! parameters to be read in from this solvers section in the sif 
  LOGICAL :: couple_to_ebfm, couple_to_icon         ! define which component is coupled to Elmer 

  CHARACTER(LEN=1024) ::  config_file, grid_dir, model_tstep
  INTEGER :: num_parts, elmer_mesh_partitions, comm_rank, comm_size, ierror
  INTEGER :: I, t, ierr
  INTEGER, POINTER :: t_icePerm(:), smbPerm(:), runoffPerm(:)
  LOGICAL :: Parallel, FirstTime=.TRUE., UnFoundFatal=.TRUE.
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Variable_t), POINTER :: t_iceVar, smbVar, runoffVar, ZsSol

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

  ! TODO: remove this check when ICON coupling is implemented
  IF (couple_to_icon) THEN
    CALL FATAL(SolverName,'>Couple To ICON< is currently not supported. Please set to FALSE')
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! retrieve the timestep in hours 
  Mesh => Solver % Mesh
  write(model_tstep, *) int(dt * 8760)
  WRITE(Message,*) 'ELMER timestep size in hours:', TRIM(model_tstep)
  CALL INFO(SolverName,Message,Level=3)
  IF (FirstTime) THEN


    ! get mesh
    ThisMesh => GetMesh(Solver)

    ! check if this is a parallel run
    IF ((ParEnv % PEs <= 1) .AND. ( .NOT. ThisMesh % SingleMesh )) THEN
      CALL FATAL(SolverName,'Only parallel runs can use this solver')
    ELSE
      grid_dir= TRIM(ThisMesh % Name)
      elmer_mesh_partitions = ParEnv % PEs
      !PRINT *, TRIM(grid_dir), elmer_mesh_partitions
      WRITE(Message,*) 'Running on ', TRIM(grid_dir),' with ',ParEnv % PEs ,' partitions' 
      CALL INFO(SolverName,Message,Level=3)
    END IF

    PRINT *, "BEFORE coupling setup"
    CALL coupling_setup(TRIM(grid_dir), elmer_mesh_partitions, TRIM(model_tstep), couple_to_ebfm, couple_to_icon)
    PRINT *, "AFTER coupling setup"


    IF (couple_to_ebfm) THEN
        ! setting up Elmer-side variables for receiving YAC variables
        ! this coul dbe replaced by an automatic picking of the names and the DOFs
        ! from the coupling-deifnitions
        ! nodal variable
        ALLOCATE(t_icePerm(Mesh % NumberOfNodes), runoffPerm(Mesh % NumberOfNodes), &
            smbPerm(GetNOFActive(Solver)))
        DO i=1,Mesh % NumberOfNodes
          t_icePerm(i) = i
          runoffPerm(i) = i
        END DO
         !This is for a element(cell) variables
        DO t=1,GetNOFActive(Solver)
          smbPerm(t) = t
        END DO
        ! nodal variables (everything that is not a flux)
        CALL DefaultVariableAdd('T_ice', dofs=1, Perm = t_icePerm)
        CALL DefaultVariableAdd('runoff', dofs=1, Perm = runoffPerm)
     
        !! element wise (cell) variable
        CALL DefaultVariableAdd('smb', dofs=1, VariableType = Variable_on_elements, Perm = smbPerm)
      
        !! get surface elevation from Elmer for the first time step
        ZsSol => VariableGet( Model % Mesh % Variables, "Zs" ,UnFoundFatal=UnFoundFatal)
        DO i=1, Mesh % NumberOfNodes
          surface_height_field(i,1) = ZsSol % Values(ZsSol % Perm(i))
        END DO
    END IF
    FirstTime = .FALSE.

    
    WRITE(Message,*) "Coupling setup with ",TRIM(grid_dir)," on ",&
         elmer_mesh_partitions, " partitions for YAC coupling done"
    CALL INFO(SolverName,Message,Level=1) 
  END IF
!!!!!!!!!! DO WE HAVE TO INITIALIZE WITH EVERY CALL ? !!!!!!!!!!!!!!
  
  IF (couple_to_ebfm) THEN
      WRITE(Message,*) 'BEFORE ELMER EBFM INTERFACE'
      CALL INFO(SolverName,Message,Level=3)
      ! couple with EBFM
      CALL elmer_ebfm_interface(ParEnv % MyPE)
      WRITE(Message,*) 'AFTER ELMER EBFM INTERFACE'
      CALL INFO(SolverName,Message,Level=3)
  !CALL elmer_icon_interface(ParEnv % MyPE)


      t_iceVar => VariableGet( Mesh % Variables, 'T_ice' )  
      smbVar => VariableGet( Mesh % Variables, 'smb' )
      runoffVar => VariableGet( Mesh % Variables, 'runoff' )
      ZsSol => VariableGet( Model % Mesh % Variables, "Zs" ,UnFoundFatal=UnFoundFatal)
      WRITE(Message,*) 'AFTER GETTING VARIABLES'
      CALL INFO(SolverName,Message,Level=3)
      IF ((.NOT.ASSOCIATED(t_iceVar)) .OR. (.NOT.ASSOCIATED(smbVar)) .OR. (.NOT.ASSOCIATED(runoffVar))) THEN
        CALL FATAL(SolverName,'Elmer variables not associated')
      END IF
      
      WRITE(Message,*) 'BEFORE WRITING NODAL VALUES'
      CALL INFO(SolverName,Message,Level=3)
       !write over values for nodes
      DO i=1, Mesh % NumberOfNodes
        t_iceVar % Values(t_iceVar % Perm(i)) = t_ice_field(i,1)
        runoffVar  % Values(runoffVar %  Perm(i)) = runoff_field(i,1)
        surface_height_field(i,1) = ZsSol % Values(ZsSol % Perm(i))
      END DO
      CALL INFO(SolverName,Message,Level=3)
      ! write over values for elements
      DO t=1, GetNOFActive(Solver)
         smbVar  % Values(smbVar %  Perm(t)) = smb_field(t,1)
      END DO
      !CALL INFO(SolverName, "Test output start", Level=1)
      !PRINT *, "Size of clt_field", SIZE(clt_field, 1),"First entry:", clt_field(1,1)
      !CALL INFO(SolverName, "Test output start", Level=1)
      !CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )
  END IF
  CALL INFO(SolverName,'Coupling step done', Level=1)
  
END SUBROUTINE YAC2Elmer

