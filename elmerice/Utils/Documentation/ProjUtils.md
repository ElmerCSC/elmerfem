# Module  ProjUtils

## General Information  
- **Module Fortran File:** ProjUtils.F90  
  
## General Description  
Module containing utility subroutines for geographical transformations, i.e. fwd (LonLat => xy) and inv. (xy => LonLat) projections.  


## Usage

In particular the subroutines **xy2LonLat(x,y,Lon,Lat)** and **LonLat2xy(Lon,Lat,x,y)**, can be used to compute longitude and latitude from x,y projectied coordinates and conversely.

To use it in your own USF or Solver:
```
<USF or SOLVER>
  USE ProjUtils

  CALL xy2LonLat(x,y,Lon,Lat)
  CALL LonLat2xy(Lon,Lat,x,y)
<>  
```

Definition of the projection and parameters in the *.sif* file, **Simulation section**

```
Simulation

! Antarctic Polar Stereographic EPSG:3031
  projection type = String "polar stereographic south"
  central_meridian = Real 0.0
  latitude_of_origin = Real -71.0

! NSIDC Sea Ice Polar Stereographic North EPSG:3014
  projection type = String "polar stereographic north"
  central_meridian = Real -45.0
  latitude_of_origin = Real 70.0

End
```

## Currently supported projections:  

- polar stereographic projections north and south, *i.e.* will cover Antarctic and Greenland applications. 
	- Implementation of the analytical solutions from J. P. Snyder, *Map-projections - A working Manual*, 1987   
	- see above for the definition of the parameters in the *.sif* file

- generic from proj4 definition:   
	- requires the fortrangis external library with proj support  
	- code is ready but currently no support for the compilation; i.e. need to update cmake

## Known Bugs and Limitations  

## Examples

An example to compute longitude and latitude using the **USF_proj.F90** user function can be found here [ELMER_TRUNK]/elmerice/Tests/Proj_South   

