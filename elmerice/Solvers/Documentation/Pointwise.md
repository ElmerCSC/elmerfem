
# Pointwise - inverse distance interolation (legacy solver!)
## General Information
- **Solver Fortran File:** pointwise.F90
- **Solver Name:** InterpolatePointValue
- **Required Output Variable(s):** main variable is a dummy - output into chosen variables from list
- **Solver Keywords:** 
  - ``Variable N`` (String) Name1 [name of variable, ``N``={1,2,3,...}]
  - ``Variable Data N`` (String) FileName1 [filename where to read xy data from]
  - ``Variable Data End N`` (String) xy [suffix where to read xy data from]
  - ``Variable DataI N`` (Integer) {0,1,2} [see meanings of switches in SIF example below]
  - ``Variable Data Offset N`` (Integer)
  - ``Variable N Supporting Points``  (Integer) 3 [minimum of data points to be used for interpolation]
  - ``Variable N Dimensions`` (Integer) 2 [dimension of variable. Here, 2, needs a file with two columns]
  - ``Variable N Exponent`` [Real] 3.0 [inverse distance weighting exponent (has to be positive) - default 2.0]
  - ``Variable N Area Scaling Factor`` (Real) 2.0 [default 1; increases the by datapoints determined max search distance by factor 2]
  - ``Variable N Directions(2)``  (Integer) 2 1 [permutation of directions for interpolation]
  
## General Description
``pointwise.F90`` was written as a way of interpolating scattered data on to the Elmer mesh based on the inverse distance method. It is included and documented due to **legacy reason, only**. Do not start using it from scratch, rather stick to [Grid2DInterpolator](./Grid2DInterpolator.md)! 


### Mulitple solutions in one SIF
Note, if pointwise is used several times within one sif-file, the corresponding .so file needs to be duplicated under a different name (see also [Flowdepth](./FlowDepth.md) Solver), i.e. ``pointwise1.so`` and ``pointwise2.so``.

## Known Bugs and Limitations
Do not use inverse distance method for extrapolation! In other words, make sure that there always is a supporting data point outside your mesh boundary.

## SIF Contents
The required keywords in the SIF file for this solver/USF are:

```
Solver 1
  Exec Solver = "Before TimeStep"
  Equation = "Pointwise Data"
  Variable =  -nooutput "Dummy"
  Variable DOFs = 1
  Procedure = "ElmerIceSolvers" "InterpolatePointValue"
  
  Variable 1 = String "mb" !Variablename in Elmer to be read in
!0-> take only Data1 as filename
!1-> "Data1"+"ceiling(Timestep*Timestepsize)+Data Offset"+Data End
!2-> "Data1"+"ceiling(Timestep*Timestepsize)+Data Offset"+"-"+ceiling(Timestep*Timestepsize)+Data Offset+1" +Data End
  Variable DataI 1 = Integer 2
  
  Variable Data 1 = String "../cmb/cmb_xyz_"  
  Variable Data Offset 1 = Integer 2005
  Variable Data End 1 = String ".dat"

  Variable 1 Supporting Points = Integer 3 !minimum of points to be used for interpolation
  Variable 1 Dimensions = Integer 2 !dimension of variable, here needs a file with two columns
  Variable 1 Exponent = Real 3.0 ! inverse distance weighting exponent (has to be positive) - default 2.0
  Variable 1 Area Scaling Factor = Real 2.0 ! increases the max search distance by factor 2 of the maximium distance of dataset poitns
  Variable 1 Directions(2) = Integer 2 1 [Here direction 1 and 2 (x and y) would be interchanged]
  
  Exported Variable 1 = -dofs 1 mb !Variablename in Elmer

End

```
## Examples
None existing




