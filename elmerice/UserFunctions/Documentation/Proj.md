# User Function Proj  
## General Information  
- **USF Fortran File:** USF_proj.F90  
- **USF Name:** xy2Lon xy2Lat LonLat2x LonLat2y  
- **Required Input Variable(s):** x,y or Lon,Lat provided as USF arguments  

## General Description  

USF_proj.F90 contains generic user functions to compute longitude and latitude from projected x,y coordinates and conversely.  
It relies on generic utilities in the module file projUtils.  

## SIF contents  
The required keywords in the SIF file for this user function are:  

e.g. to compute longitude and latitude in the initial conditions using nodal coordinates or conversely  

```
Initial Condition 1
  longitude = Variable Coordinate 1, Coordinate 2
    Real procedure "ElmerIceUSF" "xy2Lon"

  latitude = Variable Coordinate 1, Coordinate 2
    Real procedure "ElmerIceUSF" "xy2Lat"

  x = Variable longitude, latitude
    Real procedure "ElmerIceUSF" "LonLat2x"

  y = Variable longitude, latitude
    Real procedure "ElmerIceUSF" "LonLat2y"
End
```

The definition of the projectd coordinate system must be given in the **Simulation** section:
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

## Examples

