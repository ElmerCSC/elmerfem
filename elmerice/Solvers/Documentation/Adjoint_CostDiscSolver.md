## Discrete Cost function

**Module name**: Adjoint_CostDiscSolver  
**Module subroutines**: Adjoint_CostDiscSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)  
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 23.04.2020  



### Introduction
This solver computes a cost function which measures the mismatch between the model and some observations as  

If the variable is a vector we can compute the norm of the difference as (default)
$$J = \sum_{1}^{Nobs} 0.5 ||u-u^{obs}||^{2}$$  
or the difference of the norms as
$$J = \sum_{1}^{Nobs} 0.5 (||u||-||u^{obs}||)^{2}$$

If standard errors are given the cost function is scaled by the standard errors as:
$$J = \sum_{1}^{Nobs} 0.5 ||(u-u^{obs})/std||^{2}$$
or the difference of the norms as
$$J = \sum_{1}^{Nobs} 0.5 (||u||-||u^{obs}||)^{2}/||std||^2$$


Here the model variable is interpolated at the observations points and the observations do not need to be interpolated on the mesh nodes.

This solver also computes the derivative of the cost function with respect to the nodal variable. The variable that contains the derivative must exist and should be named *velocityb* if the model variable is *Flow Solution* (i.e. solving Stokes) or *ssavelocity* (i.e. solving SSA), or *Varb* for other direct equations (where *Var* is the name of the direct solver variable).

The observations can be given in an ascii or netcdf file. 

The model variable can be a scalar or a vector, and the observations can be only the first *n* components of the variable, e.g. in general only the horizontal components of the surface velocity vector are observed.

Be careful, this solver will reset the values of the cost and sensitivity to 0; so that it must be used in first place in an assimilation sequence.

In general this solver will be executed on the whole mesh for vertically integrated models or on the upper free surface
for a 3D model and 2D surface observations.

###  Keywords

Below are the related keywords in the *.sif* file:  


```
Solver *solver id* 
  
    Equation = String "Cost"  
    procedure = "ElmerIceSolvers" "Adjoint_CostDiscSolver"
    !## No need to declare a variable, this is done internally to insure 
    !## that Solver structures exist

    !# If the variable is a vector, 
    !# we can compute the difference of the norms instead the norm of the difference
     Compute norm difference = Logical [default: false]

    !# If standard errors are given they are used to normalise the cost function
     Use standard error = Logical [default: false]
    !# a min value for the standard error to avoid a division by 0
     standard error Min Value = Real [default: 100*AEPS]
     
     !# Name of an ascii file that will contain the cost value as
     !## Time, Lambda*J, RMS=sqrt(2J/(Nobs*Lambda))
     !## If Use standard error = True the 4th column will give 
     !##  the rms (i.e. not normalised by standard error)
     Cost Filename = File ""
     
     !# Name of the variable that contain the cost value
     !#  must exist and can be a global variable
     Cost Variable Name = String ""
     
     !# a multiplicatif factor can be used to scale the cost function
     Lambda=Real [default 1.0]
     
     !# The name of the model variable that is observed
     Observed Variable Name = String ""
     
     !# Name of the file that contains the observations
     !## can be netcdf (with .nc extension) or ascii (any other extension)
     !## reading netcdf is only possible for observations in a 2D plane
     !## Coordinates dimension of the Obs. must be consistent with Pb dimension;
     !##  i.e. CoordsystemDim if solver is on the bulk or 
     !###      (CoordsystemDim-1) if solver is on a boundary
     !## If the observation file is in ASCII it must containt the following columns:
     !##  coordinates(1:ndim), observations(1:vardim), standard errors(1:vardim) (optional)
     Observation File Name = File "" 
     
     !# If the variable is a vector, how many component do we observe?
     !## it will be the first *n* components
     Observed Variable dimension = Integer ...
   
     !# If true data that have been found in the mesh at the 
     !#  first visit will be saved in an ascii file
     !## will save coordinates, observation, standard error (optional), element number
     !## there will be one file per partition
     Save used data = Logical 
     !## if yes name of the output file
      output data file name = File ""
     !## and directory
      output directory = File ""
     
     
     !## Keywords related to Observations in ASCII:
     
     !# if True each partition will read her own file; partition number precede the suffix
      Parallel Observation Files = Logical  

     !# if true expect to find the Active element number in last column 
     !# (i.e. we re-read an ascii file saved by the same solver)
     Pre-Processed File = Logical 

     !# Set to True if the element number should be interpreted as the the active element number
     Element Number is Active Number = Logical 
   
     !## Keywords related to Observations in NETCDF:
     X Dim Name = File "" ![default "x"] ! name of the dimension for x
     Y Dim Name = File "" ![default "y"] ! name of the dimension for y
     
     X Var Name = File "" ![default "x"] ! name of the variable for x
     Y Var Name = File "" ![default "y"] ! name of the variable for x
     
     !# Name of the variable that contain the observation
     !## if Observed Variable dimension = 1
     Netcdf Var Name = File "" [default "Name Of The model Observed Variable"]
     !## else
     Netcdf Var 1 Name = File "" [default "vx"]
     Netcdf Var 2 Name = File "" [default "vy"]
     !# The solver will look for the attribute _FillValue to check if the data exist or not
      
     !# IF Use standard error=TRUE, the name of the variables for standard error
      Netcdf Std Var 1 Name = File "" [default "ex"]
      Netcdf Std Var 2 Name = File "" [default "ey"]
  End

```


It is possible to use a passive condition in the body force, if we want to skip the evaluation of the cost function for observations that fall within passive elements.

```
Body Force i
 !# the name of the  solver variable is NameOfEquation_var
 !# keywords relative with passive elements can be used
 Cost_var Passive = ...
End
```

### Limitations and possible improvements

- The search algorithm to locate the observations in the mesh is very efficient if the solver is executed on the whole mesh (e.g. for vertcially integrated models); however it is not efficient if the solver is executed on a boundary. 
In this case, if working with Elmer internal extrusion, it can be advantageus to :  

   - first execute this solver alone on the 2D footprint using 
      - *save used data = Logical TRUE*, 

   - then re-read the saved files using  
     - *Pre-Processed File = Logical TRUE* and  
     - *Element Number is Active Number = Logical TRUE*. 

   - If running parallel, the same number of partitions must be used and set
     - *Parallel Observation Files = Logical True*

- If the observed variable is a vector, data will only be used only if all the observed components are valid. The solver could be updated to use independently all the observed components.

Below is a list of features that are not currently possible in this solver but that could be implemented:

- For the moment we assume that errors on the observation are not correlated, i.e. the observation error covariance matrix is diagonal with the square of the standard errors in the diagonal. 

- We could allow for the possibility that only the norm of a vector is observed or only the components in a given direction.


### Tests and Examples

- See examples for the [inverse methods](../../examples/Inverse_Methods)

