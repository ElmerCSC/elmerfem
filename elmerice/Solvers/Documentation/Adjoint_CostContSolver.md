##  Continuous Cost Function

**Module name**: Adjoint_CostContSolver  
**Module subroutines**: Adjoint_CostContSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)  
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 23.04.2020  


### Introduction

This solver computes a cost function as a continuous integral over the model domain

$J = \int_{Omega} Cost$

The cost is defined in the body forces and will usually depend on some observations that need to be interpolated
at the mesh nodes.


This solver also computes the derivative of the cost function with respect to some model nodal variable that is solution of a solver  (Stokes, SSA, etc..).
The variable that contains the derivative must exist. 
To be consistent with other solvers it should be named *velocityb* if the model variable is *Flow Solution* (i.e. solving Stokes) or *ssavelocity* (i.e. solving SSA), or *Varb* for other direct equations (where *Var* is the name of the direct solver variable). 
It's dimension should be consistent with the dimension of the Model Variable. The derivative of the cost function should also be provided by the user in the Body Force and should provide the exact derivative of the Cost with respect to the given componenet of the direct variable.

Be carefull, this solver will reset the values of the cost and sensitivity to 0; so that it must be used in first place in an assimilation sequence.

In general this solver will be executed on the whole mesh for vertically integrated models or on the upper free surface
for a 3D model and 2D surface observations.


### Keywords

Bellow are the related keywords in the *.sif* file.
The solver will usually be executed on the bulk for vertically integrated models or on a boundary for 3D models.

```
Solver *solver id* 
  
    Equation = String "Cost"  
    procedure = "ElmerIceSolvers" "Adjoint_CostDiscSolver"
    ## No need to declare a variable, this is done internally to insure that Solver structures exist
     
     # Name of an ascii file that will contain the cost value as
     ## Time, J, sqrt(2J/Area)
     Cost Filename = File ""
     
     # Name of the variable that contain the cost value
     #  must exist and can be a global variable
     Cost Variable Name = String ""
     
     # The name of the model variable that contains the sensitivity
     Sensitivity Variable Name= String ""
     
      
  End

```
The cost function and its derivative must be provided in the body forces.
It is possible to use a passive condition in the body force, if we want to skip the evaluation
of the cost function in passive elements.
```
Body Force i
 Adjoint Cost = REAL ...

 # If the sensitivity Variable is a scalar:
 Adjoint Cost der = REAL ...
 # else provide the derivative with respect to each component i as:
 Adjoint Cost der i = REAL ...

 # the name of the  solver variable is NameOfEquation_var
 # keywords relative with passive elements can be used
 Cost_var Passive = ...

End
```

### Tests and Examples

