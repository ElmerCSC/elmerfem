# Test geographic coordinate system transformations

## Description

Check the cordinate transformation implemented in the ProjUtils module, in particular using the Proj library.

See documentation under elmerfem/elmerice/Utils/Documentation/ProjUtils.md 

The test :
1. read an input file with x(m),y(m),longitude(degrees),latitude(degrees)
2. Check the transformations (x,y) -> (Lon,Lat) -> (x,y)
3. Result in Fatal errors if the differences are larger than a given treshold

Currently tested projections:
- EPSG:3413 : NSIDC Sea Ice Polar Stereographic North - Using Proj and hard coded
- EPSG:3031 : Antarctic Polar Stereographic - Using Proj and hard coded
- EPSG:32632: UTM 32N
- EPSG:27572.txt: NTF (Paris) / Lambert zone II

The input files are random points in the area of interest : Mont-Blanc moutain range for EPSG:32632, EPSG:27572.txt, Greenland and Antarctica for the respective projections

Refrence Longitude and latitudes have been computed using **cs2cs** from PROJ v. 9.8.1

## Runing the test case

ElmerGrid 1 2 rectangle.grd
ElmerSolver Case.sif


## Results established:
- 02.05.2026 by Fabien Gillet-Chaulet, IGE
