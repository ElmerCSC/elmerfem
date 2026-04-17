!> Collect grid data for YAC coupling
!> Extracts grid information from mesh and converts coordinates
SUBROUTINE collect_coupling_grid_data(ThisMesh, lon_vertices, lat_vertices, &
                                      lon_cells, lat_cells, cell_to_vertex, &
                                      num_vertices_per_cell, cell_ids, vertex_ids)
  USE Types, ONLY: Mesh_t, Element_t, dp
  USE ProjUtils, ONLY: xy2LonLat, deg2rad

  IMPLICIT NONE

  TYPE(Mesh_t), POINTER :: ThisMesh
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_vertices(:), lat_vertices(:)
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_cells(:), lat_cells(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_to_vertex(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: num_vertices_per_cell(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_ids(:), vertex_ids(:)

  ! Local variables
  TYPE(Element_t), POINTER :: element
  INTEGER :: i, n, vertex_offset, v_end
  INTEGER :: nbr_vertices, nbr_cells
  INTEGER, POINTER :: this_cell_ids(:)
  
  ! Grid arrays for coupling
  REAL(KIND=dp), ALLOCATABLE :: x_vertices(:), y_vertices(:)
  REAL(KIND=dp), ALLOCATABLE :: x_cells(:), y_cells(:)

  ! Extract grid information for coupling
  nbr_vertices = ThisMesh % NumberOfNodes
  ALLOCATE(vertex_ids(nbr_vertices))
  vertex_ids = ThisMesh % ParallelInfo % GlobalDofs

  nbr_cells = ThisMesh % NumberOfBulkElements
  ALLOCATE(cell_ids(nbr_cells), num_vertices_per_cell(nbr_cells))
  DO i=1, nbr_cells
    element => ThisMesh % Elements(i)
    cell_ids(i) = element % GElementIndex
    num_vertices_per_cell(i) = element % Type % NumberOfNodes
  END DO

  ALLOCATE(cell_to_vertex(SUM(num_vertices_per_cell(:))))
  vertex_offset = 0
  DO i=1, nbr_cells
    element => ThisMesh % Elements(i)
    n = num_vertices_per_cell(i)
    v_end = vertex_offset+n
    cell_to_vertex(vertex_offset+1:v_end) = element % NodeIndexes(1:n)
    vertex_offset = v_end
  END DO

  ALLOCATE(x_vertices(nbr_vertices), y_vertices(nbr_vertices))
  x_vertices(:) = ThisMesh % Nodes % x(1:nbr_vertices)
  y_vertices(:) = ThisMesh % Nodes % y(1:nbr_vertices)

  ALLOCATE(x_cells(nbr_cells), y_cells(nbr_cells))
  DO i=1,nbr_cells
    element => ThisMesh % Elements(i)
    n = element % Type % NumberOfNodes
    this_cell_ids => element % NodeIndexes
    x_cells(i) = SUM(x_vertices(this_cell_ids(1:n))) / n
    y_cells(i) = SUM(y_vertices(this_cell_ids(1:n))) / n
  END DO

  ! Convert using ProjUtils for comparison (does NOT modify input arrays - INTENT(IN))
  ! Note: ProjUtils returns degrees, need to convert to radians for YAC
  ALLOCATE(lon_vertices(nbr_vertices), lat_vertices(nbr_vertices))
  ALLOCATE(lon_cells(nbr_cells), lat_cells(nbr_cells))
  
  DO i = 1, nbr_vertices
    CALL xy2LonLat(x_vertices(i), y_vertices(i), &
                   lon_vertices(i), lat_vertices(i))
    ! Convert from degrees to radians
    lon_vertices(i) = lon_vertices(i) * deg2rad
    lat_vertices(i) = lat_vertices(i) * deg2rad
  END DO
  
  DO i = 1, nbr_cells
    CALL xy2LonLat(x_cells(i), y_cells(i), &
                   lon_cells(i), lat_cells(i))
    ! Convert from degrees to radians
    lon_cells(i) = lon_cells(i) * deg2rad
    lat_cells(i) = lat_cells(i) * deg2rad
  END DO

  ! Clean up local arrays
  DEALLOCATE(x_vertices, y_vertices, x_cells, y_cells)

END SUBROUTINE collect_coupling_grid_data

SUBROUTINE YAC2Elmer( Model,Solver,dt,TransientSimulation )
  USE DefUtils, ONLY: GetSolverParams, GetMesh, GetNOFActive, &
    DefaultVariableAdd, GetLogical, GetLogical, GetString, ListGetString, &
    GetSimulation, MAX_NAME_LEN, VariableGet, ParEnv, variable_on_elements
  USE GeneralUtils, ONLY: I2S
  USE Types, ONLY: Model_t, Solver_t, Mesh_t, Variable_t, ValueList_t, dp, Element_t
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

  CHARACTER(LEN=1024) ::  config_file, model_tstep, coupling_timestep, grid_crs, proj_type
  INTEGER :: I, t, ierr, dt_hours, coupling_hours
  INTEGER, POINTER :: t_icePerm(:), smbPerm(:), runoffPerm(:)
  ! INTEGER, POINTER :: cltPerm(:), prPerm(:)  ! ICON is not supported at the moment
  LOGICAL :: Parallel, FirstTime=.TRUE., UnFoundFatal=.TRUE.
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Variable_t), POINTER :: t_iceVar, smbVar, runoffVar, ZsSol
  ! TYPE(Variable_t), POINTER :: cltVar, prVar  ! ICON is not supported at the moment
  REAL(KIND=dp), ALLOCATABLE :: lon_vertices(:), lat_vertices(:)
  REAL(KIND=dp), ALLOCATABLE :: lon_cells(:), lat_cells(:)
  INTEGER, ALLOCATABLE :: cell_to_vertex(:), num_vertices_per_cell(:)
  INTEGER, ALLOCATABLE :: cell_ids(:), vertex_ids(:)

  LOGICAL        :: Found

  INTERFACE
    SUBROUTINE collect_coupling_grid_data(ThisMesh, lon_vertices, lat_vertices, &
                                          lon_cells, lat_cells, cell_to_vertex, &
                                          num_vertices_per_cell, cell_ids, vertex_ids)
      USE Types, ONLY: Mesh_t, dp
      IMPLICIT NONE
      TYPE(Mesh_t), POINTER :: ThisMesh
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_vertices(:), lat_vertices(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_cells(:), lat_cells(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_to_vertex(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: num_vertices_per_cell(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_ids(:), vertex_ids(:)
    END SUBROUTINE collect_coupling_grid_data
  END INTERFACE
  
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

  ! infer grid CRS (coordinate reference system) from projection type in Simulation
  proj_type = ListGetString(GetSimulation(),'projection type',UnFoundFatal=.True.)
  SELECT CASE (TRIM(proj_type))
    CASE ('polar stereographic north')
      grid_crs = 'EPSG:3413'
    CASE ('polar stereographic south')
      grid_crs = 'EPSG:3031'
    CASE DEFAULT
      CALL FATAL(SolverName, 'Unsupported >projection type< for YAC coupling: ' // TRIM(proj_type) // &
        '. Supported values are "polar stereographic north" and "polar stereographic south".')
  END SELECT
  CALL INFO(SolverName, &
    'Using coordinate reference system (CRS): ' // TRIM(grid_crs) // ' (from projection type: ' // &
    TRIM(proj_type) // ')', Level=3)

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

    ! Collect coupling grid data (extract mesh, convert coordinates)
    CALL collect_coupling_grid_data(ThisMesh, lon_vertices, lat_vertices, &
                    lon_cells, lat_cells, cell_to_vertex, &
                    num_vertices_per_cell, cell_ids, vertex_ids)

    ! Setup YAC coupling with precomputed lon/lat coordinates (radians)
    CALL coupling_setup(lon_vertices, lat_vertices, lon_cells, lat_cells, &
              cell_to_vertex, num_vertices_per_cell, &
              cell_ids, vertex_ids, &
              TRIM(grid_crs), TRIM(ADJUSTL(I2S(coupling_hours))), &
              couple_to_ebfm, couple_to_icon)

    DEALLOCATE(lon_vertices, lat_vertices, lon_cells, lat_cells)
    DEALLOCATE(cell_to_vertex, num_vertices_per_cell, cell_ids, vertex_ids)

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

