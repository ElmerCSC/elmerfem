SUBROUTINE YAC2Elmer( Model,Solver,dt,TransientSimulation )
  USE DefUtils
  USE SolverUtils
  USE elmer_coupling
  USE elmer_icon_coupling
  
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
  INTEGER, POINTER :: cltPerm(:), prPerm(:)
  LOGICAL :: Parallel, FirstTime=.TRUE.
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Variable_t), POINTER :: cltVar, prVar

  LOGICAL        :: Found
  
  SAVE elmer_mesh_partitions, grid_dir, cltPerm, prPerm
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

    ALLOCATE(cltPerm(Mesh % NumberOfNodes),prPerm(GetNOFActive(Solver)))
    DO i=1,Mesh % NumberOfNodes
      cltPerm(i) = i
    END DO
    DO t=1,GetNOFActive(Solver)
      prPerm(t) = t
    END DO
    
    CALL DefaultVariableAdd('tas', dofs=1, Perm = cltPerm)
 
    ! element wise (cell) variable
    CALL DefaultVariableAdd('pr_snow', dofs=1, VariableType = Variable_on_elements, Perm = prPerm)
  
    
    
    FirstTime = .FALSE.

    
    WRITE(Message,*) "Coupling setup with ",TRIM(grid_dir)," on ",&
         elmer_mesh_partitions, " partitions for YAC coupling done"
    CALL INFO(SolverName,Message,Level=1) 
  END IF
!!!!!!!!!! DO WE HAVE TO INITIALIZE WITH EVERY CALL ? !!!!!!!!!!!!!!
  
  WRITE(Message,*) 'BEFORE ELMER ICON INTERFACE'
  CALL INFO(SolverName,Message,Level=3)
  ! couple with ICON
  CALL elmer_icon_interface(ParEnv % MyPE)
  WRITE(Message,*) 'AFTER ELMER ICON INTERFACE'
  CALL INFO(SolverName,Message,Level=3)


  cltVar => VariableGet( Mesh % Variables, 'tas' )  
  prVar => VariableGet( Mesh % Variables, 'pr_snow' )
  WRITE(Message,*) 'AFTER GETTING VARIABLES'
  CALL INFO(SolverName,Message,Level=3)
  IF ((.NOT.ASSOCIATED(cltVar)) .OR. (.NOT.ASSOCIATED(prVar))) THEN
    CALL FATAL(SolverName,'Elmer variables not associated')
  END IF
  
  WRITE(Message,*) 'BEFORE WRITING VALUES'
  CALL INFO(SolverName,Message,Level=3)
  ! write over values for nodes
  DO i=1, Mesh % NumberOfNodes
    !IF (ParEnv % MyPE == 0) PRINT *,i, clt_field(i,1)
    cltVar % Values(cltVar % Perm(i)) = clt_field(i,1)
    !prVar  % Values(prVar %  Perm(i)) = pr_field(i,1)
  END DO
  WRITE(Message,*) 'BEFORE WRITING SECOND VALUES'
  CALL INFO(SolverName,Message,Level=3)
  ! write over values for elements
  DO t=1, GetNOFActive(Solver)
     prVar  % Values(prVar %  Perm(t)) = pr_field(t,1)
  END DO
  !CALL INFO(SolverName, "Test output start", Level=1)
  !PRINT *, "Size of clt_field", SIZE(clt_field, 1),"First entry:", clt_field(1,1)
  !CALL INFO(SolverName, "Test output start", Level=1)
  !CALL MPI_BARRIER( ELMER_COMM_WORLD, ierr )
  CALL INFO(SolverName,'Coupling step done', Level=1)
  
END SUBROUTINE YAC2Elmer

