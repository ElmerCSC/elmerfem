# Module  ProjUtils

## General Information  
- **Module Fortran File:** ProjUtils.F90  
  
## General Description  
Module containing utility subroutines for geographical transformations, i.e. fwd (LonLat => xy) and inv. (xy => LonLat) projections.  


## Usage

In particular the subroutines **xy2LonLat(x,y,Lon,Lat)** and **LonLat2xy(Lon,Lat,x,y)**, can be used to compute longitude and latitude from x,y projected coordinates and conversely.

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

### Hard codes - Always available
- polar stereographic projections north and south, *i.e.* will cover Antarctic and Greenland applications. 
	- Implementation of the analytical solutions from J. P. Snyder, *Map-projections - A working Manual*, 1987   
	- see above for the definition of the parameters in the *.sif* file

### Using the PROJ library

- can be any projection supported by the PROJ library:
    - see https://proj.org/en/stable/
	- requires to compile Elmer with PROJ (see below)
    - Uses the PROJ 6 API

- .sif content:   

```
Simulation

  projection type = String "proj"
  ! the current mesh CRS definition
  !  can be every format supported by proj. ESPG codes are the preferred definitions
  CRS code = File "EPSG:32632"

End
```

- Compilation using CMAKE:  

```
...
 -DWITH_PROJ:BOOL=TRUE \
 -DPROJ_LIBRARIES="${PATH_TO_PROJ_LIB_DIR}/libproj.so" \
...
```

[!WARNING]  
 A previous version was relying on the FortranGIS interface and the PROJ 4 API. As of June 2026, this is no more supported 
  and the configuration and compilation must be done as described above

## Known Bugs and Limitations  

## Examples and Tests

- An example to compute longitude and latitude using the **USF_proj.F90** user function can be found here [ELMER_TRUNK]/elmerice/Tests/Proj_South   

- Various usual projections are tested under [ELMER_TRUNK]/elmerice/Tests/Proj6:  
    - EPSG:3031 - Antarctic : Hard Coded and using Proj
    - EPSG:3413 - Greenland : Hard Coded and using Proj
    - EPSG:32632 - UTM 32N : using Proj
    - EPSG:27572 - Lambert II : using Proj

