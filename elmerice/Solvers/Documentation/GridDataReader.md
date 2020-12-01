#GridDataReader Solver
##General Information
- **Solver Fortran File:** GridDataReader.f90 (under elmerice/netcdf2-tree)
- **Solver Name:** GridDataReader
- **Required Output Variable(s):** None
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** all from NetCDF file interpolated variables
- **Optional Input Variable(s):** None

##Compilation
**Update 37f46d0** under elmerice branch
Now you can automatically build, test & compile the solver with cmake.
**Update e375a11** under elmerice branch:
cmake variable WITH_GridDataReader no more required. GridDatareader it compiled if netcdf is found.

KeyWords:
```
cmake
   -DNetCDF_INCLUDE_DIR=PATH_TO_INCLUDE \
   -DNetCDF_LIBRARY=PATH_TO_libnetcdf.so \
   -DNetCDFF_LIBRARY=PATH_TO_libnetcdff.so \
```
###Old Version
This solver is not part of the standard Elmer/Ice installation and has to be compiled manually. There is an example compilation call in the directory. Naturally, the system needs to have a compatible NetCDF library installed. In particular, the Fortran compiler with which the NetCDF library has been compiled has to be the same Elmer and Elmer/Ice has been compiled with (module compatibility)

##General Information
This auxiliary solver enables to read in variables from a NetCDF file and interpolates the values to the mesh as variables of the same name. The interpolation is done by using the Finite Element test functions.

By default the code will loop over the mesh nodes and read only the values required for the interpolation at the given node. This might involve a lot of I/O and become a bottle neck. Loading the whole netcdf at the begining requires to allocate more space in the memory but might be more efficient if possible. This possibility is given by setting: *Read full array = Logical True*

##Known Bugs
None, so far.

##SIF contents
```
Solver 1
  Equation = Reader
  Procedure = "GridDataReader" "GridDataReader"

  !---- NOTE: File is case sensitive, String is not!
  Filename = File "../greenland.nc" 
  
  !----- Load the whole netcdf array instead of inquiring only the values required 
  !  for the interpolation at the current mesh node
  Read full array = Logical [default: false]

  Time Dim Name = String "time"
  X Dim Name = String "X" 
  Y Dim Name = String "Y"

  Time Var Name = String "time"
  X Var Name = String "x" 
  Y Var Name = String "y"

  !--- Interpolation variables tolerances
  X Epsilon = Real 1.0e-2 
  Y Epsilon = Real 1.0e-2 
  Time Epsilon = Real 0.01

  Interpolation Bias = Real 0.0 
  Interpolation Multiplier = Real 1.0 

  Is Time Counter = Logical True

  Variable 1 = annualtemp
  Variable 2 = julytemp
  Variable 3 = preciptation
  Variable 4 = ablation
  
  Enable Scaling = Logical True ! Scales the Elmer grid to match the NetCDF grid - dangerous
End
```

##Test
Located under [ELMER_TRUNK]/elmerice/Tests/GridDataReader
You can test with

```
ctest -L GridDataReader
or 
ctest -L netcdf
```

##Example
The directory [ELMER_TRUNK]/elmerice/netcdf2 contains the example above, that should work with the NetCDF file provided by the [searise group](http://websrv.cs.umt.edu/isis/index.php/Present_Day_Greenland#Data_Download).
