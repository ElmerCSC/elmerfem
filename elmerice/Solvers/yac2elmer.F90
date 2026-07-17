!> Collect grid data for YAC coupling
!> Extracts grid information from mesh and converts coordinates
SUBROUTINE collect_coupling_grid_data(ThisMesh, lon_vertices, lat_vertices, &
                                      lon_cells, lat_cells, cell_to_vertex, &
                                      num_vertices_per_cell, cell_ids, &
                                      vertex_ids, boundary_corner_mask)
  USE Types, ONLY: Mesh_t, Element_t, dp
  USE ProjUtils, ONLY: xy2LonLat, deg2rad
  USE DefUtils, ONLY: GetBoundaryEdgeIndex
  USE Messages, ONLY: FATAL
  USE MeshUtils, ONLY: FindMeshEdges

  IMPLICIT NONE

  TYPE(Mesh_t), POINTER :: ThisMesh
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_vertices(:), lat_vertices(:)
  REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_cells(:), lat_cells(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_to_vertex(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: num_vertices_per_cell(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_ids(:), vertex_ids(:)
  LOGICAL, ALLOCATABLE, INTENT(OUT) :: boundary_corner_mask(:)

  ! Local variables
  TYPE(Element_t), POINTER :: element
  INTEGER :: i, n, vertex_offset, v_end, bnd_elem_idx, boundary_edge_idx
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
  ALLOCATE(boundary_corner_mask(nbr_vertices))
  boundary_corner_mask = .FALSE.

  IF (ThisMesh % MeshDim /= 2) THEN
    CALL FATAL('collect_coupling_grid_data', &
      'boundary_corner_mask is implemented for 2D meshes only')
  END IF

  DO i=1, nbr_cells
    element => ThisMesh % Elements(i)
    cell_ids(i) = element % GElementIndex
    num_vertices_per_cell(i) = element % Type % NumberOfNodes
  END DO

  ! Ensure edge tables are built (needed for GetBoundaryEdgeIndex)
  IF (.NOT. ASSOCIATED(ThisMesh % Edges)) CALL FindMeshEdges(ThisMesh)

  ! Mark only true physical domain-boundary cells; partition-only boundaries stay false.
  DO bnd_elem_idx = ThisMesh % NumberOfBulkElements + 1, &
                    ThisMesh % NumberOfBulkElements + ThisMesh % NumberOfBoundaryElements
    element => ThisMesh % Elements(bnd_elem_idx)
    ! Ensure that all boundary elements are edges (should be the case for 2D meshes)
    IF (element % Type % NumberOfEdges /= 1) THEN
      CALL FATAL('collect_coupling_grid_data', &
        'Expected exactly one edge per 2D boundary element')
    END IF
    ! Filter to actual boundary edges only.
    IF (.NOT. ASSOCIATED(element % BoundaryInfo)) CYCLE

    boundary_edge_idx = GetBoundaryEdgeIndex(element, 1)
    IF (boundary_edge_idx <= 0) THEN
      CALL FATAL('collect_coupling_grid_data', &
        'GetBoundaryEdgeIndex returned invalid edge index for boundary element')
    END IF

    ! Check if this is a parallel mesh
    IF (ASSOCIATED(ThisMesh % ParallelInfo % EdgeInterface)) THEN
      IF (boundary_edge_idx > SIZE(ThisMesh % ParallelInfo % EdgeInterface)) THEN
        CALL FATAL('collect_coupling_grid_data', &
          'Boundary edge index exceeds ParallelInfo % EdgeInterface size')
      END IF
      ! check if boundary edge is a partition boundary edge
      IF (.NOT. ThisMesh % ParallelInfo % EdgeInterface(boundary_edge_idx)) THEN
        boundary_corner_mask(element % NodeIndexes(1:element % Type % NumberOfNodes)) = .TRUE.
      END IF
    ELSE
      ! For a non-parallel mesh all boundary edges are physical boundaries
      boundary_corner_mask(element % NodeIndexes(1:element % Type % NumberOfNodes)) = .TRUE.
    END IF
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
                   lon_vertices(i), lat_vertices(i), .False.)
    ! Convert from degrees to radians
    lon_vertices(i) = lon_vertices(i)
    lat_vertices(i) = lat_vertices(i)
  END DO

  DO i = 1, nbr_cells
    CALL xy2LonLat(x_cells(i), y_cells(i), &
                   lon_cells(i), lat_cells(i), .False.)
    ! Convert from degrees to radians
    lon_cells(i) = lon_cells(i)
    lat_cells(i) = lat_cells(i)
  END DO

  ! Clean up local arrays
  DEALLOCATE(x_vertices, y_vertices, x_cells, y_cells)

END SUBROUTINE collect_coupling_grid_data

SUBROUTINE YAC2Elmer( Model,Solver,dt,TransientSimulation )
  USE DefUtils, ONLY: GetSolverParams, GetMesh, GetNOFActive, &
    DefaultVariableAdd, GetLogical, GetString, GetConstReal, &
    ListGetString, ListGetConstReal, &
    GetSimulation, MAX_NAME_LEN, VariableGet, ParEnv, variable_on_elements, &
    GetActiveElement
  USE GeneralUtils, ONLY: I2S
  USE Types, ONLY: Model_t, Solver_t, Mesh_t, Variable_t, ValueList_t, dp, Element_t
  USE Messages, ONLY: Message, FATAL, INFO, USE_YAC
  USE elmer_coupling, ONLY: coupling_setup, is_root_rank
  USE elmer_ebfm_coupling, ONLY: elmer_ebfm_interface, t_ice_field, smb_field, &
                                 runoff_field, surface_height_field
  USE elmer_icon_coupling, ONLY: elmer_icon_interface, t_oce_post_field, &
                                 sal_oce_post_field, liquid_ice_sheet_flux_field
  USE MPI

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

  CHARACTER(LEN=1024) ::  config_file, grid_crs, proj_type
  CHARACTER(LEN=1024) ::  coupling_timestep
  REAL(KIND=dp) :: coupling_timestep_in_years, local_flux_sum, global_flux_sum
  CHARACTER(LEN=1024) ::  yac_calendar, yac_start_time, yac_end_time
  INTEGER :: i, t, ierr
  INTEGER, POINTER :: t_icePerm(:), smbPerm(:), runoffPerm(:)
  INTEGER, POINTER :: t_ocePerm(:), sal_ocePerm(:)
  REAL(KIND=dp) :: central_meridian, latitude_of_origin
  REAL(KIND=dp) :: expected_central_meridian, expected_latitude_of_origin
  CHARACTER(LEN=16) :: expected_central_meridian_str, expected_latitude_of_origin_str
  LOGICAL :: Parallel, FirstTime=.TRUE., UnFoundFatal=.TRUE.
  TYPE(Mesh_t),POINTER :: Mesh
  TYPE(Variable_t), POINTER :: t_iceVar, smbVar, runoffVar, ZsSol
  TYPE(Variable_t), POINTER :: t_oceVar, sal_oceVar, bmb_fluxVar
  TYPE(Element_t), POINTER :: Element
  REAL(KIND=dp), ALLOCATABLE :: lon_vertices(:), lat_vertices(:)
  REAL(KIND=dp), ALLOCATABLE :: lon_cells(:), lat_cells(:)
  INTEGER, ALLOCATABLE :: cell_to_vertex(:), num_vertices_per_cell(:)
  INTEGER, ALLOCATABLE :: cell_ids(:), vertex_ids(:)
  LOGICAL, ALLOCATABLE :: boundary_corner_mask(:)

  ! Variables needed for user output
  CHARACTER(LEN=1024) :: coupling_timestep_in_years_str
  CHARACTER(LEN=1024) :: dt_str

  LOGICAL        :: Found

  INTERFACE
    SUBROUTINE collect_coupling_grid_data(ThisMesh, lon_vertices, lat_vertices, &
                                          lon_cells, lat_cells, cell_to_vertex, &
                                          num_vertices_per_cell, cell_ids, vertex_ids, &
                                          boundary_corner_mask)
      USE Types, ONLY: Mesh_t, dp
      IMPLICIT NONE
      TYPE(Mesh_t), POINTER :: ThisMesh
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_vertices(:), lat_vertices(:)
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT) :: lon_cells(:), lat_cells(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_to_vertex(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: num_vertices_per_cell(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: cell_ids(:), vertex_ids(:)
      LOGICAL, ALLOCATABLE, INTENT(OUT) :: boundary_corner_mask(:)
    END SUBROUTINE collect_coupling_grid_data
  END INTERFACE

  SolverParams => GetSolverParams()

  ! read config file
  config_file = GetString(SolverParams, 'Config File Name',  Found )
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Config File Name< found in yac2elmer solver" &
    )
  END IF

  ! read YAC calendar
  yac_calendar = GetString(SolverParams, 'Calendar', Found)
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Calendar< found in yac2elmer solver" &
    )
  END IF
  yac_calendar = TRIM(ADJUSTL(yac_calendar))
  SELECT CASE (yac_calendar)
    CASE ("proleptic_gregorian", "year_360", "year_365")
      CONTINUE
    CASE DEFAULT
      CALL FATAL(SolverName, &
        "Unsupported >Calendar< value: '" // TRIM(yac_calendar) // "'. " // &
        "Accepted values are proleptic_gregorian, year_360, year_365." &
      )
  END SELECT

  ! read YAC coupling start time (ISO 8601)
  yac_start_time = GetString(SolverParams, 'Coupling Start Time', Found)
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Coupling Start Time< found in yac2elmer solver" &
    )
  END IF
  IF (LEN_TRIM(yac_start_time) == 0) THEN
    CALL FATAL(SolverName, &
      "Keyword >Coupling Start Time< must not be empty" &
    )
  END IF

  ! read YAC coupling end time (ISO 8601)
  yac_end_time = GetString(SolverParams, 'Coupling End Time', Found)
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Coupling End Time< found in yac2elmer solver" &
    )
  END IF
  IF (LEN_TRIM(yac_end_time) == 0) THEN
    CALL FATAL(SolverName, &
      "Keyword >Coupling End Time< must not be empty" &
    )
  END IF

  ! check if config file actually exists
  INQUIRE(FILE=TRIM(config_file), EXIST=USE_YAC)
  IF (.NOT. USE_YAC) THEN
    CALL FATAL(SolverName, &
      "Coupling config file not found" &
    )
  END IF

  ! check if it is coupled to EBFM
  couple_to_ebfm = GetLogical(SolverParams, 'Couple To EBFM',  Found )
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Couple To EBFM< found in yac2elmer solver" &
    )
  END IF

  ! check if it is coupled to ICON
  couple_to_icon = GetLogical(SolverParams, 'Couple To ICON',  Found )
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Couple To ICON< found in yac2elmer solver. " &
    )
  END IF

  ! read coupling timestep as ISO 8601 string
  coupling_timestep = GetString(SolverParams, &
    "Coupling Time Step", Found)
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Coupling Time Step< found in yac2elmer solver. " // &
      "IMPORTANT: User has to ensure that >Coupling Time Step< and " // &
      ">Coupling Time Step in yrs< are consistent." &
    )
  END IF
  IF (LEN_TRIM(coupling_timestep) == 0) THEN
    CALL FATAL(SolverName, &
      "Keyword >Coupling Time Step< must not be empty. " // &
      "IMPORTANT: User has to ensure that >Coupling Time Step< and " // &
      ">Coupling Time Step in yrs< are consistent." &
    )
  END IF

  ! read coupling timestep in years (used for Elmer-internal consistency checks)
  !
  ! This is a duplication of 'Coupling Time Step' and the user is responsible to
  ! ensure that both values are consistent. The reason for this duplication is
  ! that the coupling timestep in years is needed for Elmer-internal consistency
  ! checks, while the ISO 8601 string is needed for YAC coupling setup.
  !
  ! Conversion of ISO 8601 string to years is not trivial and would require
  ! additional parsing logic that is able to consider edge cases (leap years).

  coupling_timestep_in_years = GetConstReal(SolverParams, &
    "Coupling Time Step in yrs", Found)
  IF (.NOT. Found) THEN
    CALL FATAL(SolverName, &
      "No keyword >Coupling Time Step in yrs< found in yac2elmer solver. " // &
      "IMPORTANT: User has to ensure that >Coupling Time Step< and " // &
      ">Coupling Time Step in yrs< are consistent." &
    )
  END IF
  IF (LEN_TRIM(coupling_timestep) == 0) THEN
    CALL FATAL(SolverName, &
      "Keyword >Coupling Time Step in yrs< must not be empty. " // &
      "IMPORTANT: User has to ensure that >Coupling Time Step< and " // &
      ">Coupling Time Step in yrs< are consistent." &
    )
  END IF

  ! infer grid CRS (coordinate reference system) from projection type in Simulation
  proj_type = ListGetString(GetSimulation(),'projection type',UnFoundFatal=.True.)

  SELECT CASE (TRIM(proj_type))
    CASE ('polar stereographic north')
      grid_crs = 'EPSG:3413'
      expected_central_meridian = -45.0_dp
      expected_latitude_of_origin = 70.0_dp
    CASE ('polar stereographic south')
      grid_crs = 'EPSG:3031'
      expected_central_meridian = 0.0_dp
      expected_latitude_of_origin = -71.0_dp
    CASE DEFAULT
      CALL FATAL(SolverName, 'Unsupported >projection type< for YAC coupling: ' // TRIM(proj_type) // &
        '. Supported values are "polar stereographic north" and "polar stereographic south".')
  END SELECT

  ! Do consistency checks
  central_meridian = ListGetConstReal(GetSimulation(), 'central_meridian', UnFoundFatal=.True.)
  latitude_of_origin = ListGetConstReal(GetSimulation(), 'latitude_of_origin', UnFoundFatal=.True.)

  ! Format expected values as strings
  WRITE(expected_central_meridian_str, '(F6.1)') expected_central_meridian
  WRITE(expected_latitude_of_origin_str, '(F6.1)') expected_latitude_of_origin

  IF (ABS(central_meridian - expected_central_meridian) > 1.0e-6_dp) THEN
    CALL FATAL(SolverName, 'For >projection type< "' // TRIM(proj_type) &
      // '", >central_meridian< must be ' // TRIM(expected_central_meridian_str))
  END IF
  IF (ABS(latitude_of_origin - expected_latitude_of_origin) > 1.0e-6_dp) THEN
    CALL FATAL(SolverName, 'For >projection type< "' // TRIM(proj_type) &
      // '", >latitude_of_origin< must be ' // TRIM(expected_latitude_of_origin_str))
  END IF

  CALL INFO(SolverName, &
    "Using coordinate reference system (CRS): " // TRIM(grid_crs) // " " // &
    "(from projection type: " // TRIM(proj_type) // ")", Level=3 &
  )

  IF (.NOT. (couple_to_ebfm .OR. couple_to_icon)) THEN
    CALL FATAL(SolverName, &
      "At least one of >Couple To EBFM< or >Couple To ICON< must be TRUE" &
    )
  END IF

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Mesh => Solver % Mesh

  WRITE(coupling_timestep_in_years_str,'(G0)') coupling_timestep_in_years
  WRITE(dt_str,'(G0)') dt

  ! Check if coupling_timestep_in_years is greater than or equal to dt
  IF (coupling_timestep_in_years < dt) THEN
    CALL FATAL( SolverName, &
      "'Coupling Time Step' must be greater than or equal to Elmer time step size " // &
      "('Coupling Time Step in yrs': " // TRIM(coupling_timestep_in_years_str) // &
      " and Elmer dt = " // TRIM(dt_str) // ")" &
    )
  END IF
  ! Check if coupling_timestep_in_years is an integer multiple of dt
  IF (MOD(coupling_timestep_in_years, dt) > 1.0e-6_dp) THEN
    CALL FATAL( SolverName, &
      "'Coupling Time Step' must be an integer multiple of Elmer time step size " // &
      "('Coupling Time Step in yrs': " // TRIM(coupling_timestep_in_years_str) // &
      " and Elmer dt = " // TRIM(dt_str) // ")" &
    )
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
                    num_vertices_per_cell, cell_ids, vertex_ids, &
                    boundary_corner_mask)

    ! Setup YAC coupling with precomputed lon/lat coordinates (radians)
    CALL coupling_setup(lon_vertices, lat_vertices, lon_cells, lat_cells, &
              cell_to_vertex, num_vertices_per_cell, &
              cell_ids, vertex_ids, &
              TRIM(grid_crs), &
              TRIM(config_file), &
              TRIM(yac_calendar), &
              TRIM(yac_start_time), &
              TRIM(yac_end_time), &
              TRIM(coupling_timestep), &
              couple_to_ebfm, couple_to_icon, &
              boundary_corner_mask)

    DEALLOCATE(lon_vertices, lat_vertices, lon_cells, lat_cells)
    DEALLOCATE(cell_to_vertex, num_vertices_per_cell, cell_ids, vertex_ids)
    DEALLOCATE(boundary_corner_mask)

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
      ALLOCATE(t_ocePerm(Mesh % NumberOfNodes), sal_ocePerm(Mesh % NumberOfNodes))
      DO i=1,Mesh % NumberOfNodes
        t_ocePerm(i) = i
        sal_ocePerm(i) = i
      END DO
      CALL DefaultVariableAdd('temp_oce', dofs=1, Perm = t_ocePerm)
      CALL DefaultVariableAdd('sal_oce', dofs=1, Perm = sal_ocePerm)
      ! Initialize bmb_flux_field for first time step
      bmb_fluxVar => VariableGet( Model % Mesh % Variables, "bmb_flux", UnFoundFatal=UnFoundFatal)
      DO t=1, GetNOFActive(Solver)
        liquid_ice_sheet_flux_field(t,1) = bmb_fluxVar % Values(bmb_fluxVar % Perm(t))
      END DO
    END IF

    FirstTime = .FALSE.

    CALL INFO(SolverName, "YAC coupling setup done", Level=1)
  END IF
!!!!!!!!!! DO WE HAVE TO INITIALIZE WITH EVERY CALL ? !!!!!!!!!!!!!!
  IF (couple_to_icon) THEN
      CALL INFO(SolverName, 'Starting coupling Elmer -> ICON-O', Level=3)
      CALL INFO(SolverName, 'Getting Elmer liquid and solid flux variables', Level=30)
      ! Update bmb_flux_field before sending to ICON
      bmb_fluxVar => VariableGet( Model % Mesh % Variables, "bmb_flux", UnFoundFatal=UnFoundFatal)
      DO t=1, GetNOFActive(Solver)
        liquid_ice_sheet_flux_field(t,1) = bmb_fluxVar % Values(bmb_fluxVar % Perm(t))
      END DO
      ! Write total liquid ice sheet flux to log
      local_flux_sum = 0.0_dp
      DO t = 1, GetNOFActive(Solver)
        Element => GetActiveElement(t, Solver)
        IF (ParEnv % myPe /= Element % partIndex) CYCLE
        local_flux_sum = local_flux_sum + liquid_ice_sheet_flux_field(t,1)
      END DO
      CALL MPI_Allreduce(local_flux_sum, global_flux_sum, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, ParEnv % ActiveComm, ierr)
      WRITE(Message,'(A,F15.3)') 'Sum of liquid_ice_sheet_flux_field: ', global_flux_sum
      CALL INFO(SolverName, Message, Level=3)
      ! couple with ICON-O
      CALL elmer_icon_interface(is_root_rank)
      CALL INFO(SolverName, 'Finished coupling Elmer -> ICON-O', Level=3)
      CALL INFO(SolverName, 'Writing variables from ICON-O into Elmer variables...', Level=30)
      t_oceVar => VariableGet( Mesh % Variables, 'temp_oce' )
      sal_oceVar => VariableGet( Mesh % Variables, 'sal_oce' )
      IF ((.NOT.ASSOCIATED(t_oceVar)) .OR. (.NOT.ASSOCIATED(sal_oceVar))) THEN
        CALL FATAL(SolverName,'Elmer variables not associated')
      END IF
      ! write over values for nodes
      DO i=1, Mesh % NumberOfNodes
        t_oceVar % Values(t_oceVar % Perm(i)) = t_oce_post_field(i,1)
        sal_oceVar % Values(sal_oceVar % Perm(i)) = sal_oce_post_field(i,1)
      END DO
  END IF

  IF (couple_to_ebfm) THEN
      CALL INFO(SolverName, 'Starting coupling Elmer -> EBFM', Level=3)
      ! couple with EBFM
      CALL elmer_ebfm_interface(is_root_rank)
      CALL INFO(SolverName, 'Finished coupling Elmer -> EBFM', Level=3)
      CALL INFO(SolverName, 'Writing variables from EBFM into Elmer variables...', Level=30)
      t_iceVar => VariableGet( Mesh % Variables, 'T_ice' )
      smbVar => VariableGet( Mesh % Variables, 'smb' )
      runoffVar => VariableGet( Mesh % Variables, 'runoff' )
      ZsSol => VariableGet( Model % Mesh % Variables, "Zs" ,UnFoundFatal=UnFoundFatal)
      IF ((.NOT.ASSOCIATED(t_iceVar)) .OR. (.NOT.ASSOCIATED(smbVar)) .OR. (.NOT.ASSOCIATED(runoffVar))) THEN
        CALL FATAL(SolverName,'Elmer variables not associated')
      END IF

      ! Validate element variable permutations. In a parallel run, element
      ! variables declared via "Exported Variable = -elem" in the SIF get a
      ! proper parallel Perm (built by SetActiveElementsTable). Without that
      ! declaration, the variable may end up with an invalid Perm.

      ! ! "smb" is known to require this declaration
      IF (MINVAL(smbVar % Perm(1:GetNOFActive(Solver))) <= 0) THEN
        CALL FATAL(SolverName, &
          'smb variable has zero/invalid Perm entries. ' // &
          'Add  Exported Variable = -elem "smb"  to the YAC2Elmer solver block in the SIF.')
      END IF

      ! "runoff" did not cause issues, but adding the same check for safety
      IF (MINVAL(runoffVar % Perm(1:GetNOFActive(Solver))) <= 0) THEN
        CALL FATAL(SolverName, &
          'runoff variable has zero/invalid Perm entries. ' // &
          'Add  Exported Variable = -elem "runoff"  to the YAC2Elmer solver block in the SIF.')
      END IF

      CALL INFO(SolverName, 'Writing nodal variables...', Level=30)
       !write over values for nodes
      DO i=1, Mesh % NumberOfNodes
        t_iceVar % Values(t_iceVar % Perm(i)) = t_ice_field(i,1)
        surface_height_field(i,1) = ZsSol % Values(ZsSol % Perm(i))
      END DO
      CALL INFO(SolverName, 'Writing cell variables....', Level=30)
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

  CALL INFO(SolverName,'Coupling step done', Level=1)

END SUBROUTINE YAC2Elmer

