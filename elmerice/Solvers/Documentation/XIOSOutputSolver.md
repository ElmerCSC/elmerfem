# Solver XIOSOutputSolver
## General Information
- **Solver  Fortran File:** XIOSOutputSolver.F90   
- **Solver/User Function Name:** XIOSOutputSolver  
- **Required Output Variable(s):** None  
- **Required Input Variable(s):** None  
- **Optional Output Variable(s):** None  
- **Optional Input Variable(s):** None
- **Solver Keywords:** 
  - Scalar Field i : Name of nodal and elemntal variables that can be requested by XIOS, i=1,...  
  - Global Variable i: Name of a global variable that can be requested by XIOS, e.g. time
  
## General Description

Interface to [XIOS](https://forge.ipsl.jussieu.fr/ioserver), a library designed to manage NETCDF outputs of climate models.

XIOS supports writing unstructured data using the [UGRID convention](http://ugrid-conventions.github.io/ugrid-conventions/#ugrid-conventions-v10). This is currently restricted to 2D.

### Getting XIOS

1. Get XIOS, it requires the **trunk** version and it is not working with the latest (2.5) release. 
```
svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/trunk XIOS
```  

2. Compile XIOS, see its documentation   
   - You will have to add **-fPIC** for the *BASE_CFLAGS* and *BASE_FFLAGS* as we wil need a shared library
   - To use the UGRID format we need to have netcdf with parallel I/O support


3. Create the shared library; XIOS compilation only creates a static library **libxios.a** but we require a shared library. 
locate the share library then do:   
```
# ! extract .o files
ar -x libxios.a  
#! make a shared library from the .o; update with the c-compilator you have been using, e.g. mpiicc for intel
mpicc -shared *.o -o libxios.so
#! remove the .o
rm -f *.o
```   
The directory where you store the shared libary should be in your *LD_LIBRABRY_PATH* environnement variable, if in a non standard location.

### Compiling Elmer with XIOS   

If you have XIOS, in general you will also compile ELMER with Netcdf support, so the corresponding libraries should already by provided with the corresponding CMAKE instructions.    
**Compile Elmer with the same NETCDF and compilor as XIOS!**  


Regarding NETCDF and HDF5, you should used shared libraries:
Assuming NETCDF is installed under *NETCDF_DIR* and HDF5 under *HDF5_DIR*
Adding the following instructions to your cmake configuration file should be sufficient.


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

Assuming XIOS is installed under *XIOS_DIR*:
Adding the following instructions to your cmake configuration file should be sufficient.

```
...
 -DWITH_XIOS:BOOL=TRUE \
 -DXIOS_INCLUDE_DIR=${XIOS_DIR}/include -DXIOS_ROOT=${XIOS_DIR}  \
 -DXIOS_LIBRARY="${XIOS_DIR}/lib/libxios.so" \
...
```   


### Running Elmer with XIOS

Elmer will call XIOS if Elmer has been compiled with XIOS support and the XIOS configuration file **iodef.xml** is present in the current directory.  The context *id* should be **elmerice** to configure the outputs.

It can be run in attached mode, where each Elmer MPI process will call one instance of xios, using
```
mpirun -np N ElmerSolver_mpi
```
or in detached mode using
```
mpirun -np N ElmerSolver_mpi : -np N2 xios_server
```

## Known Bugs and Limitations

- Restricted to 2D, but should be possible to extand this to a 3D simulation if working with a vertically extruded mesh, and we want to save only 2D variables.

- Should work with 303 or 404 elements, *to check for mesh with a mixture how to prescribe the bounds?*

- Should work with *halo elements* as they are skipped for the saving table

- should work with higher order elements, e.g. 306, as we will save only the corners. *To check how to get the max number of corners*

- XIOS automatically re-computes the connectivity tables, including edges. 
	- *Can we directly provide this to XIOS?*
	- *To see how to save variable on edges, e.g. grounding line flux..."*

## SIF Contents
The required keywords in the SIF file for this solver are

```
Simulation
..
!# Definition of the projection, see ProjUtils documentation
! Greenland EPSG:3413
  projection type = String "polar stereographic north"
  central_meridian = Real -45.0
  latitude_of_origin = Real 70.0
End


Solver ..
  ! Usually executed after time steps
   Exec Solver = After Timestep

   Equation = "XIOSOutPutSolve"
   Procedure = "ElmerIceSolvers" "XIOSOutputSolver"

! node and elem vars
   Scalar Field i = String "VarName" ! Variable/component Name

  !Global Variables
   Global Variable i = String "Global Variable Name"

End

```

## elmerice context

For the XIOS xml configuration file, the context *id* should be **elmerice**.   
A variable with **id="time_units** should be provided to define the units of the time step, i.e. **1y** if we are using years or **1s** for seconds. The time step send to xios is then the Elmer time step *dt* times the *time_units*.    

- :warning: **its is to the user responsability to check that the time step is constant and a finite fraction of the output frequency**   
- :warning: **id for the variable in the xml file should correspond to the Elmer variable name provided in the .sif file but in lower case**, i.e. in the .sif file "VarName" is case isensitive and should be referred as **id="varname"** for XIOS.

```
<context id="elmerice">

<!-- define the time unit system... should be lenght 2 -->
<variable_definition>
	<variable id="time_units" type="string">1y</variable>
</variable_definition>

<!-- id corespond to var names as defined in the .sif file and should be lower case!! -->
<field id="varname"  name=... />

</context>
```


## Visualising the resulting UGRID file.	

UGRID Netcdf files can be visualized with:   
- QGIS as **mesh layers** using the [crayfish plugin](https://www.lutraconsulting.co.uk/projects/crayfish)   
	- the dimension and associated variable **time_counter** must be changed to **time**:
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
	- As the coordinates are in lon/lat, it might be interesting to provide the projected coorinates and thus change the names of the coordinate in the mesh attributes? or improve the reader?

- [psyplot](https://psyplot.github.io/)
	- not adapted yet for nodal variable as it re-compute a delaunay triangulation?

- [gridded](https://github.com/NOAA-ORR-ERD/gridded) 
	- not tested

## Examples
An example can be found here [ELMER_TRUNK]/elmerice/Tests/XIOS

