## Regularisation solver {#adjoint_CostReg}

**Module name**: Adjoint_CostRegSolver  
**Module subroutines**: Adjoint_CostRegSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)  
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 24.04.2020  


### Introduction
This solver computes a cost function classically used for regularisation of the inverse problems 
and it's derivative with respect to the input nodal variable **V**.

Two regularisations are possible:  

- First option penalises 1st spatial derivatives *[default]*:

> $$ J_{reg} = \int_{\Omega} 0.5  (|dV/dx|)^2 d\Omega $$

The dimension of the spatial derivatives is *CoordinateSystemDimension()* or *CoordinateSystemDimension()-1*
depending if the solver is executed on the bulk on a boundary.

- Second option penalises difference from a *prior*:

> $$ J_{reg} = \int_{\Omega} 0.5  ((V-V^{prior})/s^2)^2 d\Omega $$


where $V^{prior}$ is the *prior* estimate and $s^2$ the variance 


### Keywords

Bellow are the related keywords in the *.sif* file:  


```
Solver *solver id* 
  
    Equation = String "CostReg"  
    procedure = "ElmerIceSolvers" "Adjoint_CostRegSolver"
    !## No need to declare a variable, this is done internally to insure
    !    that Solver structures exist
 
     !# True if cost function and gradient must be initialised 
     !   to 0 in this solver
     Reset Cost Value = Logical [Default: True]

     !# Name of an ascii file that will contain the cost value as
     !## Time, J_{reg}, RMS=sqrt(2J_{reg}/Area)
     Cost Filename = File ""
     
     !# Name of the variable that contain the cost value
     !#  must exist and can be a global variable
     Cost Variable Name = String ""
     
     !# a multiplicatif factor can be used to scale the cost function
     Lambda=Real [default 1.0]
     
     !# The name of the model variable that will contain the derivative
     !  of J w.r.t. the input variable
     Gradient Variable Name = String ""

     !# The name of the model variable V
     !## if not found will look for the keyword "CostReg Nodal Variable" 
     !    in the body force.
     Optimized Variable Name = String ""

     !# Switch to regularisation from *prior*
     A priori Regularisation = Logical [default:False]
     
      
  End

```

Values for the nodal variable **V** can be given in the body force if  *Optimized Variable Name* was
not provided in the solver parameters. 
This can be usefull if we want to apply some change of variable compared to the variable that is optimised,
However it will be to the user responsability to provide the exact derivative of his change of variable 
as, by default, the solver computes the derivative with respect to the nodal value.
In this cas ethe user can directly provide the derivative with the keyword *CostReg Nodal Variable der = Real ...*
see below and *[examples]*

If using *a priori* regularisation, nodal values for the *prior* and *standard deviation* are given in the body force too.

```
Body Force i
 # value for the nodal *V*
  CostReg Nodal Variable = Real ...
 # If the definition above implies a change of variable
 #  the derivative of the function above can be performed using 
 #  the following keyword
  CostReg Nodal Variable derivative = Real ...

 # value for the nodal *V^{prior}*
  CostReg Nodal Prior = Real ...
 # value for the nodal std *s*
  CostReg Nodal std = Real ...
End
```

It is possible to use a passive condition in the body force, if we want to skip the evaluation of the cost function within passive elements.

```
Body Force i
 # the name of the  solver variable is NameOfEquation_var
 # keywords relative with passive elements can be used
 CostReg_var Passive = ...
End
```

### Limitations and possible improvments

Bellow is a list of features that are not currently possible in this solver but that could be implemented:  

- If running on a boundary, we use only the first spatial derivatives, so this solver is intended to be used on a bottom or top surface.

- If Regularisation with respect to the *prior* is used, for the moment we assume no spatial error in the statitics and use only a standard deviation; The cost function could easily be improved if a full background error covariance matrix is known


### Tests and Examples

Validation examples available under:  

- ELMER_SRC/elmerice/examples/Adjoint_CostRegSolver

   Tested:

   - Regularisation penalising first spatial derivatives on a 2D mesh and on the bottom boundary of a 3D mesh 

   - Regularisation from  *prior*

   - Optimisation  

   - Prescribing nodal value from the body forces using a change of variable


- See examples for the [inverse methods](../../examples/Inverse_Methods)
