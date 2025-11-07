SUBROUTINE YAC2Elmer( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE SolverUtils
  USE elmer_coupling
  USE elmer_ebfm_coupling
  
  IMPLICIT NONE

  TYPE(Model_t)  :: Model
  TYPE(Solver_t) :: Solver
  REAL(KIND=dp)  :: dt
  LOGICAL        :: TransientSimulation

  
  TYPE(ValueList_t), POINTER :: SolverParams
  TYPE(Mesh_t), POINTER :: ThisMesh
  CHARACTER(LEN=MAX_NAME_LEN):: SolverName='YAC2Elmer'
  CHARACTER(LEN=1024) ::  config_file, grid_dir, model_tstep
  INTEGER :: num_parts, elmer_mesh_partitions, comm_rank, comm_size, ierror
  INTEGER :: I, t, ierr
  INTEGER, POINTER :: t_icePerm(:), smbPerm(:), runoffPerm(:)
  LOGICAL :: Parallel, FirstTime=.TRUE.
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Variable_t), POINTER :: t_iceVar, smbVar, runoffVar

  LOGICAL        :: Found
  
  SAVE elmer_mesh_partitions, grid_dir, t_icePerm, smbPerm, runoffPerm
  ! check if config file exist. Terminate program if it is not found!

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

  Mesh => Solver % Mesh
  write(model_tstep, *) int(dt * 8760)
  WRITE(Message,*) 'ELMER timestep size in hours:', TRIM(model_tstep)
  CALL INFO(SolverName,Message,Level=3)
  
  !!!!!!!!!! DO WE HAVE TO INITIALIZE WITH EVERY CALL ? !!!!!!!!!!!!!!
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

    !CALL coupling_setup(TRIM(grid_dir), elmer_mesh_partitions, "1")
    PRINT *, "BEFORE coupling setup"
    CALL coupling_setup(TRIM(grid_dir), elmer_mesh_partitions, TRIM(model_tstep))
    PRINT *, "AFTER coupling setup"


    ! setting up Elmer-side variables for receiving YAC variables
    ! this coul dbe replaced by an automatic picking of the names and the DOFs
    ! from the coupling-deifnitions
    ! nodal variable

    ALLOCATE(t_icePerm(Mesh % NumberOfNodes), runoffPerm(Mesh % NumberOfNodes), smbPerm(Mesh % NumberOfNodes))
    DO i=1,Mesh % NumberOfNodes
      t_icePerm(i) = i
      runoffPerm(i) = i
      smbPerm(i) = t
    END DO
    ! This would be for a element(cell) variables
    !DO t=1,GetNOFActive(Solver)
      !smbPerm(t) = t
    !END DO

    PRINT *, "AFTER allocation"
    
    ! nodal variables (everything that is not a flux)
    CALL DefaultVariableAdd('T_ice', dofs=1, Perm = t_icePerm)
    CALL DefaultVariableAdd('smb', dofs=1, Perm = smbPerm)
    CALL DefaultVariableAdd('runoff', dofs=1, Perm = runoffPerm)
 
    !! element wise (cell) variable
    !CALL DefaultVariableAdd('smb', dofs=1, VariableType = Variable_on_elements, Perm = smbPerm)
    !CALL DefaultVariableAdd('runoff', dofs=1, VariableType = Variable_on_elements, Perm = runoffPerm)
  
    PRINT *, "AFTER variable addition"
    
    
    FirstTime = .FALSE.

    
    WRITE(Message,*) "Coupling setup with ",TRIM(grid_dir)," on ",&
         elmer_mesh_partitions, " partitions for YAC coupling done"
    CALL INFO(SolverName,Message,Level=1) 
  END IF
!!!!!!!!!! DO WE HAVE TO INITIALIZE WITH EVERY CALL ? !!!!!!!!!!!!!!
  
  WRITE(Message,*) 'BEFORE ELMER EBFM INTERFACE'
  CALL INFO(SolverName,Message,Level=3)
  ! couple with EBFM
  CALL elmer_ebfm_interface(ParEnv % MyPE)
  !CALL elmer_icon_interface(ParEnv % MyPE)
  WRITE(Message,*) 'AFTER ELMER EBFM INTERFACE'
  CALL INFO(SolverName,Message,Level=3)


  t_iceVar => VariableGet( Mesh % Variables, 'T_ice' )  
  smbVar => VariableGet( Mesh % Variables, 'smb' )
  runoffVar => VariableGet( Mesh % Variables, 'runoff' )
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
    smbVar  % Values(smbVar %  Perm(t)) = smb_field(t,1)
    runoffVar  % Values(runoffVar %  Perm(i)) = runoff_field(i,1)
  END DO
  WRITE(Message,*) 'BEFORE WRITING CELL VALUES'
  CALL INFO(SolverName,Message,Level=3)
  ! write over values for elements
  !DO t=1, GetNOFActive(Solver)
     !smbVar  % Values(smbVar %  Perm(t)) = smb_field(t,1)
     !!runoffVar  % Values(runoffVar %  Perm(t)) = runoff_field(t,1)
  !END DO
  !CALL INFO(SolverName, "Test output start", Level=1)
  !PRINT *, "Size of clt_field", SIZE(clt_field, 1),"First entry:", clt_field(1,1)
  !CALL INFO(SolverName, "Test output start", Level=1)
  !CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )
  CALL INFO(SolverName,'Coupling step done', Level=1)
  
END SUBROUTINE YAC2Elmer

