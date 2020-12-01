#Grid2DInterpolator
##General Description
This solver interpolates data given on a regular 2D grid in an ASCII file (x y Value). A bilinear interpolation is used. By default, the data in the ASCII file have to be ordered such that

```
  x1 y1 val11
  x2 y1 val21
  ...
  xn y1 valn1
  x1 y2 val12
  ...
  xn ym valnm 
```
If the second column (the y-value) is changing faster, i.e.,

```
  x1 y1 val11
  x1 y2 val12
  ...
  x1 ym val1m 
  ...
  xn y1 valn1
  xn y2 valn2
  ...
  xn ym valnm
```
one can give the keyword

 `Variable 1 Invert = Logical True`
In any case, it is essential, that the order (x1 < x2 < x3 < … < xn) as well as (y1 < y2 < y3 < … < ym) applies.

Further, if there are those points with no-data simply missing from the file (saves space), the keyword

 `Variable 1 Fill = Logical True`
will insert either the default no-data value or the one given with the Variable 1 No Data value at the missing coordinates. Mind, that if this keyword is either False or not given, any missing entry will stop the routine.

The grid is described by giving:
- (x0, y0) the left-bottom corner coordinate
- (lx, ly) the x and y lengths of the covered domain
- (Nx, Ny) the number of levels in x and y directions (n,m in the tables above)
- Cells with no data are identified by value (within a given tolerance)
WARNING: All data with values between (noData - tol) and (noData + tol) will be ignored (where noData and tol are the no data value and tolerance, see below). In the case a node falls in a cell of the DEM containing no data, the nodal value is assigned to the closest DEM value. This robust feature might be dangerous if too many nodes are outside of the DEM domain. The solver gives the number of nodes outside of the DEM domain and the farthest DEM point used to evaluate the nodal values.

##SIF contents
Add this solver and execute it before the simulation. Here, it is used to read the surface DEM and the bedrock DEM given by the files DEM_TR_surf.dat and DEM_TR_bed.dat, respectively.

```
Solver 1
  Exec Solver = Before Simulation
  Equation = "Read DEM"
  Procedure = "ElmerIceSolvers" "Grid2DInterpolator"

  Variable 1 = String "ZsDEM"
  Variable 1 data file = File "DEM_TR_surf.dat"
  Variable 1 x0 = REal 947700.0d0
  Variable 1 y0 = REal 2104850.0d0
  Variable 1 lx = REal 800.0
  Variable 1 ly = REal 350.0
  Variable 1 Nx = Integer 268
  Variable 1 Ny = Integer 118
  Variable 1 no data = Real -999.0
  Variable 1 no data tol = Real 0.1

  Variable 2 = String "bedrockDEM"
  Variable 2 data file = File "DEM_TR_bed.dat"
  Variable 2 x0 = REal 947700.0d0
  Variable 2 y0 = REal 2104850.0d0
  Variable 2 lx = REal 600.0
  Variable 2 ly = REal 350.0
  Variable 2 Nx = Integer 301
  Variable 2 Ny = Integer 176
End
```
Note that the “no data” and “no data tol” keywords are optional for each variable. If these keywords are not given in the sif file the no data value defaults to -9999.0 and the tolerance to 0.001.

The variables bedrockDEM and ZsDEM have to exported from an other solver

```
Solver 3
  Equation = "Navier-Stokes"
  Exported Variable 1 = -dofs 1 "ZsDEM"
  Exported Variable 2 = -dofs 1 "bedrockDEM"
  ...
End
```
The variables bedrockDEM and ZsDEM can then be used to estimate the Bottom Surface and Top Surface keywords used by the StructuredMeshMapper solver (have a look [here](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:structuredmesh)). Here, a minimum ice thickness of 1 meter is imposed.
```
! Bedrock BC
Boundary Condition 2
  Bottom Surface = Equals bedrockDEM
  ...
End

! Upper Surface
Boundary Condition 3
  Top Surface = Variable ZsDEM, bedrockDEM
    Real MATC "if (tx(0)>tx(1)+1.0) {tx(0)} else {tx(1)+1.0}"
End
```

##Example
An example using the Tete Rousse surface and bedrock DEM can be found in [ELMER_TRUNK]/elmerice/Tests/Grid2DInterpolator.
