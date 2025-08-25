# Solver XIOSOutputSolver
## General Information
- **Solver  Fortran File:** XIOSOutputSolver.F90
- **Solver/User Function Name:** XIOSOutputSolver
- **Required Output Variable(s):** None
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None
- **Solver Keywords:**
  - Scalar Field i : Name of nodal and elemental variables that can be requested by XIOS, i=1,...
  - Global Variable i: Name of a global variable that can be requested by XIOS, e.g. time

## General Description

Interface to [XIOS](https://forge.ipsl.fr/ioserver), a library designed to manage NETCDF outputs of climate models.

XIOS supports writing unstructured data using the [UGRID convention](http://ugrid-conventions.github.io/ugrid-conventions/#ugrid-conventions-v10). This is currently restricted to 2D.

### Getting XIOS

1. Get the **trunk** version of XIOS (Elmer does not work with the latest 2.5 release):
```
svn co https://forge.ipsl.fr/ioserver/svn/XIOS/trunk XIOS
```

2. Compile XIOS (please refer to its documentation), but note that:
   - You will have to add **-fPIC** for the *BASE_CFLAGS* and *BASE_FFLAGS* as Elmer needs a shared library.
   - To use the UGRID format we need to have netcdf with parallel I/O support.

3. Create the shared XIOS library. By default, the XIOS compilation process only creates a static library **libxios.a**, but Elmer requires a shared library. Locate the static library then do:
```
# Extract the object files (*.o) from the static library
ar -x libxios.a
# Create the shared library from the object files (use the same compiler you used to compile XIOS!)
mpicc -shared *.o -o libxios.so
# Remove the object files
rm -f *.o
```
The directory which contains the shared library should be listed in your *LD_LIBRARY_PATH* environment variable if it is not a standard location.

### Compiling Elmer with XIOS

If you have XIOS, in general you will also compile Elmer with Netcdf support, so the corresponding libraries should already by provided with the corresponding CMAKE instructions.

:warning: **Compile Elmer with the same NETCDF and compiler as XIOS!**

You also need shared versions of the NetCDF and HDF5 libraries. Assuming NetCDF is installed under *NETCDF_DIR* and HDF5 under *HDF5_DIR*, adding the following instructions to your cmake configuration file should work:

```
...
 -DWITH_NETCDF:BOOL=TRUE \
 -DNETCDF_INCLUDE_DIR=${NETCDF_DIR}/include -DNETCDF_ROOT=${NETCDF_DIR}  \
 -DNETCDF_LIBRARY="${NETCDF_DIR}/lib/libnetcdf.so" \
 -DNETCDFF_LIBRARY="${NETCDF_DIR}/lib/libnetcdff.so"\
 -DPHDF5_INCLUDE_DIR=${HDF5_DIR}/include -DPHDF5_ROOT=${HDF5_DIR}  \
 -DPHDF5_LIBRARY="${HDF5_DIR}/lib/libhdf5.so" \
 -DPHDF5HL_LIBRARY="${HDF5_DIR}/lib/libhdf5_hl.so" \
...

```

Assuming XIOS is installed under *XIOS_DIR*, adding the following instructions to your cmake configuration file should work:

```
...
 -DWITH_XIOS:BOOL=TRUE \
 -DXIOS_INCLUDE_DIR=${XIOS_DIR}/include -DXIOS_ROOT=${XIOS_DIR}  \
 -DXIOS_LIBRARY="${XIOS_DIR}/lib/libxios.so" \
...
```

### Running Elmer with XIOS

Elmer will call XIOS if Elmer has been compiled with XIOS support and the XIOS configuration file **iodef.xml** is present in the work directory at run-time. The context *id* in the XIOS configuration file(s) should be properly set (more details below).

XIOS can be run in attached mode, where each Elmer MPI process will call one instance of XIOS, using
```
mpirun -np N ElmerSolver_mpi
```
or in detached mode using
```
mpirun -np N ElmerSolver_mpi : -np N2 xios_server
```

### Saving variables with XIOS

- To use the UGRID format, geographical coordinates are required (Longitude/Latitude). They are computed using generic functionalities from [ProjUtils](../../Utils/Documentation/ProjUtils.md). The projection description should be provided in the *Simulation* section (only south and north polar stereographic projections are supported at the moment). The projected coordinates can be saved by requesting to save the Elmer coordinates.

- Elmer variables that have to be transferred to XIOS are defined with the Solver keywords *Scalar Field i* for nodal and elemental variables and *Global Variable i* for global variables. The *id* in the XIOS configuration file should be the Elmer variable name :warning: **in lower case**.

- Nodal variables can be avaraged by element using the keyword *Scalar Field i compute cell average = Logical True*. In this case the corresponding *id* for XIOS should be *varname_elem* (:warning: **in lower case**). Multiplying by the element area and summing over the elements will provide a conservative alternative to integrating nodal variables with Elmer.

- The following variables are computed by the solver and can be accessed by XIOS:

  * The elements' surface area can be accessed with the *id* *cell_area*.

## SIF Contents

The required keywords in the SIF file for this solver are:

```
Simulation
..
 ! Definition of the projection, see ProjUtils documentation
 ! Greenland EPSG:3413
   projection type = String "polar stereographic north"
   central_meridian = Real -45.0
   latitude_of_origin = Real 70.0

 ! setting the time_step to an integer number of days (15) for a 360 day calendar.
   Output  Intervals(1) = 24
   Timestep Intervals(1) = 240
   Timestep Sizes(1) = $ 15.0/360.0
End


Solver ..
 ! Usually executed after time steps
   Exec Solver = After Timestep

   Equation = "XIOSOutPutSolve"
   Procedure = "ElmerIceSolvers" "XIOSOutputSolver"

 ! keywords related to calendar management
 ! time_units (mandatory): set the time unit system used by Elmer
   time_units=String "1y"

 ! time-step (optional): the duration of the time step; other time_step=time_units*dt
   timestep=String "15d"
 ! for consistency, check that taking 1/dt time_step leads
 ! to the same duration than time_units with xx seconds
   timestep tolerance = Real 1.0

 ! to set the start date from Elmer; start date will be reference date + (Gettime()-dt)*time_units
 ! i.e. for restart if time=10+dt => start_date=2025-01-01
   reference date=String "2015-01-01"

 ! automatically add this suffix to all files and file_groups that start by filei, i=1,9 or 01,09, and 10,99
   file names suffix = String "_$name$_$suffix$"

 ! node and elem vars
   Scalar Field i = String "VarName" ! Variable/component Name
 ! to compute cell avaraged values
   Scalar Field i compute cell average = Logical True

 ! Global Variables
   Global Variable i = String "Global Variable Name"

 ! (for debbuging) can decrease the info level (default 4) will show all requested variables by xios...
   Solver info level = integer 3
End

```

## XIOS context and execution ID

In the XIOS xml configuration file(s), the context *id* should match the value of the command-line argument `-exec-id` passed when starting the Elmer solver (its default value is *elmerice* when `-exec-id` is omitted). For example, if starting Elmer with `ElmerSolver_mpi greenland.sif -exec-id greenland`, then the context *id* should be *greenland* in XIOS' configuration file(s).

- :warning: To run multiple Elmer simulations in a single mpirun command using a single detached XIOS instance, each Elmer simulation should be given a different `-exec-id`.
- :warning: According to XIOS' documentation: **Value of unit may be integer or floating (not recommended)**, so setting the time_step to dt*"1y" might lead to spurious effect. It is then possible to provide the time-setp in the solver as a duration; We will check that a full year gives the same duration in seconds within a given tolerance as in general an integer number of days will not lead to a finite fraction of year.... Might be better to move to units in days instead of years....?
- :warning: **it is the user's responsibility to check that the time step is constant and a finite fraction of the output frequency**.
- :warning: **id for the variable in the xml file should correspond to the Elmer variable name provided in the .sif file but in lower case**, i.e. in the .sif file "VarName" is case insensitive and should be referred as **id="varname"** for XIOS. There is a sanity check that a variable defined with the keywords *Global Variable* and *Scalar Field* are defined in the XIOS configuration file.

```
<!-- mandatory context definition (use default value or the value of -exec-id) -->
<context id="elmerice">

<!-- id correspond to var names as defined in the .sif file and should be lower case!! -->
<field id="varname"  name=... />

<!-- if you want to compute element-averaged values from nodal values; add "_elem" to the var name -->
<field id="varname_elem"  name=... />

<!-- setting elmer/ice version as a file global attribute -->
<file id=... >
  <field ... />
  <!-- global attribute definition ... All variables with this id will be updated to contain both rev. and vers. numbers -->
  <variable id="elmerversion" name="model_version" type="string"> elmer ice v9.0</variable>
</file>

<!-- mandatory domains and grids. NB: the domain ids are hard-coded in Elmer but the domain name can be changed (it must be the same for the three domains though) -->
<domain_definition>
  <!-- mandatory domains ...  -->
  <domain id="cells" name="mesh2D"/>
  <domain id="edges" name="mesh2D"/>
  <domain id="nodes" name="mesh2D"/>
  <!-- ...  -->
</domain_definition>

<grid_definition>
  <!-- mandatory grids... -->
  <grid id="GridCells">
     <domain domain_ref="cells"/>
  </grid>

  <grid id="GridNodes">
    <domain domain_ref="nodes"/>
  </grid>

  <grid id="GridEdges">
    <domain domain_ref="edges"/>
  </grid>
  <!-- ...  -->
</grid_definition>

</context>
```

### Hard coded definitions

A *field_group* with **id="mesh_info"** will be automatically added is not already present, with the attribute *operation="once"*.
If not already defined the following fields (related to mesh information) will be added to the group:
- node_x: node x coordinate
- node_y: node y coordinate
- cell_area: element area
- boundary_condition: Edge % BoundaryInfo % Constraint

i.e. this is equivalent to having the following definitions in the elmerice context; and the field_group and fields can be used to be saved in files.

:warning: these variables are not recomputed; however they might change in some applications... especially the boundary condition if we have passive/active boundary conditions...

```
<field_group id="mesh_info"  operation="once" >
   <field id="node_x"  name="x"  standard_name="projection_x_coordinate"  unit="m" grid_ref="GridNodes" />
   <field id="node_y"  name="y"  standard_name="projection_y_coordinate"  unit="m" grid_ref="GridNodes" />
   <field id="cell_area" name="cell_area" unit="m2" grid_ref="GridCells" />
   <!-- boundary condition number for edges... should be better to output as int. but qgis do not support variable that are not float? -->
   <field id="boundary_condition" name="boundary_condition" unit="1" default_value="0" prec="4" grid_ref="GridEdges"/>
</field_group>
```

## Reading an unstructured Netcdf File

Files produced with XIOS can be read with the [UGridDataReader](UGridDataReader.md)

## Visualising the resulting UGRID file.

UGRID Netcdf files can be visualized with:
- QGIS as **mesh layers** using the [crayfish plugin](https://www.lutraconsulting.co.uk/projects/crayfish)
	- If the time dimension and associated variable is **time_counter**,this must be changed to **time**:
	```
        ncrename -d time_counter,time -v time_counter,time output_ugrid.nc
	````

- Paraview using the [UGRID Reader](https://github.com/FeliciaBrisc/UGRID-Reader-for-ParaView)
	- Must have a compiled version of paraview; restricted to 5.6 for the moment?
	- Actually the reader will not understand variables with no time dimension and will crash because of some missing attributes?
	We can extract the mesh topology and required variables:
	```
	ncks -v mesh2D,mesh2D_face_nodes,mesh2D_node_x,mesh2D_node_y,MyNodalVar,MyElemVar output_ugrid.nc tmp.nc
	```
	- As the coordinates are in lon/lat, it might be interesting to provide the projected coordinates and thus change the names of the coordinate in the mesh attributes? or improve the reader?

- [psyplot](https://psyplot.github.io/)
	- not adapted yet for nodal variable as it re-computes a delaunay triangulation?

- [gridded](https://github.com/NOAA-ORR-ERD/gridded)
	- not tested

### Tips

- Variables defined with the *id="elmerversion"* will be automatically updated to *Elmer/Ice vVERSION_NUMBER (Rev: REVISION_NUMBER)*, so that it can be added as a global attribute in output files.

- By default the files will contain a time dimension named **time_counter**, and the associated variable is the **time_centered** variable. As most software, e.g. for visualisation, will look for a dimension named **time**, the default can be changed using the following keywords in the file definition: *time_counter_name="time" time_counter="instant"*. But remember, the true time for a variable can be *time_instant* or *time_centered* depending on the time operator; this is defined in the variable attribute.

## Known Bugs and Limitations

- Restricted to applications where you can define geographic coordinates.

- Restricted to 2D, but should be possible to extend this to a 3D simulation if working with a vertically extruded mesh, and we want to save only 2D variables.

- Should work with 303 or 404 elements, *to check for mesh with a mixture how to prescribe the bounds?*

- Should work with *halo elements* as they are skipped from the saving table.

- should work with higher order elements, e.g. 306, as we will save only the corners?
	- *To check how to get the max number of corners?*
	- Might not be as simple as we rely on the global DOFs for the ordering...

- XIOS automatically re-computes the connectivity tables, including edges.
	- *Can we directly provide this to XIOS?*

- Edges: *Element="n:0 e:1"* is automatically added to the solver parameters so that Elmer computes the edge table and ordering. :warning: Variables expoted in this solver will be defined by edges... To see how to define edge variable in other solvers (e.g. grounding line flux) to send them to Xios...

## Examples

An example can be found here [ELMER_TRUNK]/elmerice/Tests/Xios

An example to compute element averaged-values and do the reduction with XIOS can be found in [ELMER_TRUNK]/elmerice/Tests/Xios2

An example to read variables stored in a file produced with XIOS can be found in [ELMER_TRUNK]/elmerice/Tests/UGridDataReader
