# UGRID Data Reader 
## General Information
- **Solver Fortran File:** UGridDataReader.f90 
- **Solver Name:** UGridDataReader
- **Required Output Variable(s):** None
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** variables read from the input netcdf file
- **Optional Input Variable(s):** None

## General Information

This solver reads variables in an unstructured netcdf file (e.g. following the UGRID format) at node and element locations.
It can be used to *e.g.* read variables that have been produced with the **XIOSOutPutSolver** or that have been conservatively interpolated on the mesh using e.g. cdo.

The input file structure should correspond to the **current serial mesh**, and variables should be arranged using the node and element ordering.

## Time Index

For variable that have a time dimension, the time index to read can be defined using the following keywords:
- *Is Time Counter = Logical True [Default: False]*: the time index start at 1 and is incremented by 1 at each visit. So if we have e.g. annual values taht we want to read from the begining and update every year, we can execute this solver every 1/dt time steps.

- *Time Point = Real ...*: Define a Time point to read. The Time index is then defined as *TimePoint = floor(time-dt/2) + 1*, i.e. Index=1 for Time Point in ]0,1]; Index=2 for Time Point in ]1,2], etc...

- **By default** the current simulation time is used as time point, i.e. equivalent to:
```
Time Point = Variable Time
  REAL MATC "tx[0]"
```

## Known Bugs and Limitations

- Netcdf should have the same structure as the mesh; might be interesting to add a functionnality to read a 2D file for an vertically extruded 3D mesh?

- In parallel the same file will be accessed by all partitions. Might become a bottleneck for very large simulations?

## SIF contents
```
 Solver 1
   Equation = "UGridDataReader"
   Procedure = "ElmerIceSolvers" "UGridDataReader"

   ! If the file contain a time dimesion
   ! the time step to read can be defined using the following keywords:

   !  ! time counter: increment time by 1 at each visit
   Is Time Counter = Logical [Default: False]
   !  ! set a time point
   Time Point = Real ...
   ! default use the simulation time equivalent to 
   Time Point = Real Time
     REAL MATC "tx[0]"

   ! File Name (case sensitive)
   File Name = File "output_ugrid.nc"

   ! List of variables to read
   ! Variable name as in the Netcdf (case sensitive)
   Variable Name 1 = File "MyNodalVar"
   ! the Elmer variable (must exist, as a solver or exported variable)
   Target Variable 1 = string "NodeVar"
  
  ...
End
```

## Compilation

You need to have ELMER compiled with **NETCDF**, eventualy compiled with the **HDF5** libraries. 
Assuming NETCDF is installed under *NETCDF_DIR* and HDF5 under *HDF5_DIR*
Adding the following instructions to your cmake configuration file should be sufficient.

```
cmake
...
 -DWITH_NETCDF:BOOL=TRUE \
 -DNETCDF_INCLUDE_DIR=${NETCDF_DIR}/include \
 -DNETCDF_LIBRARY="${NETCDF_DIR}/lib/libnetcdf.so" \
 -DNETCDFF_LIBRARY="${NETCDF_DIR}/lib/libnetcdff.so"\
 -DPHDF5_INCLUDE_DIR=${HDF5_DIR}/include  \
 -DPHDF5_LIBRARY="${HDF5_DIR}/lib/libhdf5.so" \
 -DPHDF5HL_LIBRARY="${HDF5_DIR}/lib/libhdf5_hl.so" \
...

```

## Test
A validation test is located under  [ELMER_TRUNK]/elmerice/Tests/UGridDataReader


