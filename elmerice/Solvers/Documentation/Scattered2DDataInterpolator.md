# Scattered2DDataInterpolator

## General Information
- **Solver Fortran File:** Scattered2DDataInterpolator.F90
- **Solver Name:** Scattered2DDataInterpolator
- **Required Output Variable(s):** None
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None
- **Solver Keywords:** see *SIF Contents*
  
## General Description

Solver to interpolate sparse 2D data sets onto the FE mesh.

Interpolation is made using external libraries  
- [nn](https://github.com/sakov/nn-c): a c code for Natural Neighbours interpolation  
- [csa](https://github.com/sakov/csa-c) : a c code for cubic spline approximation  

The user is invited to get familiar with these libraries by using the command line utilities and test case provided with the libraries

Data are read from an **ASCII file** (x,y,variable) or from a **netcdf file** (if data file has extension .nc)

As this solver depends on several external libraries (nn, csa, netcdf) it is not compiled and included by default in the ElmerIceSolver shared library.

**To use the solver :**  you have to compile the required external libraries and compile elmer with the following cmake arguments

```
cmake ...
   -DWITH_ScatteredDataInterpolator:BOOL=TRUE \
   -DNetCDF_INCLUDE_DIR=PATH_TO_INCLUDE \
   -DNetCDF_LIBRARY=PATH_TO_libnetcdf.so \
   -DNetCDFF_LIBRARY=PATH_TO_libnetcdff.so \
   -DCSA_LIBRARY=PATH_TO_libcsa.a \
   -DCSA_INCLUDE_DIR=PATH_TO_csa \
   -DNN_INCLUDE_DIR=/PATH_TO_nn \
   -DNN_LIBRARY=PATH_TO_libnn.a
```

## Known Bugs and Limitations
- None

## SIF contents
The Solver section must include:

```
Solver 
 ! This solver is intended to be used once before simulation 
 ! to import data sets onto the FE mesh  
  Exec Solver = Before simulation

  Equation = "ScatteredInter"
  procedure = "Scattered2DDataInterpolator" "Scattered2DDataInterpolator"
  

!!  Bounding Box
 ! If this parameter is set only the data points that are 
 ! within Max/Min mesh corrdinates + the real Value
 ! are used. It an be useful in parallel if all the data are stored in one file
  Bounding Box dx = Real 1.0e3  
 
!! CheckNaN 
 ! check if interpolation method returns NaN 
  CheckNaN = Logical True ! [Default: True]; 
 ! By default replace NaN by nearest available value
 ! We can replace NaNs by Real value with following flag
  Replace NaN by = Real -9999999 

!!!!! NNI or linear (nn-c library)

 ! Default Sibson interpolation
  Variable 1 = String "ZsNNI"
  Variable 1 data file = File "Rand200.txt"

!! Valid Min/Max values
!!  valid min/max values can be set with the following keywords
!!    in this case data which are below (above) the value are set to the Min (Max) value
!!    and the results will also be limited by these values
!!  Only the min or the max can be set.
   Variable 1  Valid Min Value = Real ...
   Variable 1  Valid Max Value = Real ...

  Variable 2 = String "ZsNNIW"
  Variable 2 data file = File "Rand200.txt"
  Variable 2 W = Real 0.
! W restricts extrapolation by assigning minimal allowed
!    weight for a vertex (normally "-1" or so; lower
!    values correspond to lower reliability; "0" means
!    no extrapolation)
! Default W=-HUGE(RealNumber);i.e. extrapolation allowed

  Variable 3 = String "ZsNNINS"
  Variable 3 data file = File "Rand200.txt"
  Variable 3 method = String "ns"
 ! method Non-Sibsonian interpolation (nn-c); W can be applied here too

  Variable 4 = String "ZsLin"
  Variable 4 data file = File "Rand200.txt"
  Variable 4 method = String "li"
 ! method linear interpolation (nn-c); W no effect here

  Variable 5 = String "ZsFlight"
  Variable 5 data file = File "FlightLines.txt"

!!!!! Cubic spline (csa-c library)
  Variable 6 = String "ZsCS"
  Variable 6 data file = File "Rand200.txt"
  Variable 6 method = String "cs"
 ! method cubic spline (csa)
  Variable 6 nppc = integer 5
!set the average number of points per cell (default = 5,
!works best for uniform data. Decrease to get smaller
!               cells or increase to get larger cells)
  Variable 6 k = integer 140
! set the spline sensitivity (default = 140, reduce to get
!                     smoother results)

End
```
Keywords specific to netcdf:

```
Solver 2
  Exec Solver = Before simulation
  Equation = "ScatteredInter"
  procedure = "Scattered2DDataInterpolator" "Scattered2DDataInterpolator"

 ! name of the variable in the netcdf file as it is case sensitive use "File"
  Variable 1 = File "NetcdfVar"  
! results can be stroed in an Elmer variable with a different name
! set the name of the Elmer variable with the following flag
 Target Variable 1 = String "MyVar"

 ! .nc means that the data file is in netcdf format
  Variable 1 data file = File "<NETCDF_FILE>.nc"  
  
  ! FillValue; 
  ! The solver reads the attribute _FillValue and exclude the corresponding data.
  ! It can be changed (or set if the attribute is not present) with the following flag;
  Variable 1 Fill Value = Real ... 
  
  !! Dim and Var names
  Variable 1 x-dim name = File "x" !name of the x dimension (default: x)
  Variable 1 y-dim name = File "y" !name of the y dimension (default: y)
  Variable 1 x-Var name = File "x" !name of the variable with x-coordinates (default: x)
  Variable 1 y-Var name = File "y" !name of the variable with x-coordinates (default: y)
  
  !!! keywords for the interpolation are identical. see above
  
End
```

## Test
Located under [ELMER_TRUNK]/elmerice/Tests/ScatteredDataInterpolator
You can test with

ctest -L ScatteredDataInterpolator
or 
ctest -L netcdf

## Example
A test case is provided under
[ELMER_TRUNK]/elmerice/examples/Test_Scattered2DDataInterpolator .
