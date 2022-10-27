## AdjointThickness : Gradient Solver {#thickness_gradient_solver}

**Module name**: AdjointThickness_GradientSolver.F90    
**Module subroutines**: AdjointThickness_GradientSolver     
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)      
**Document authors**: Fabien Gillet-Chaulet      
**Document edited**: 10/12/2020    

**Required input variables:**

 - Variable solution of the direct problem ([Thickness solver](#thickness_direct_solver))
 - The adjoint variable computed from the [adjoint linear system](#adjoint_linearsolver)

**output variables: [Optional]**

 - Nodal derivatives with respect to:

    - convection velocity $u$: [DJDuv]
    - the top surface accumulation $Ms$: [DJDsmbTop]
    - the bottom surface accumulation $Mb$: [DJDsmbBot]


### Introduction

This is the adjoint code of the [Thickness solver](#thickness_direct_solver).

Computes the derivatives of a cost function that involves the thickness $H$, solution of the direct solver, with respect to the input parameters (the convection velocity and the mass balance forcing).

The sequence in the .sif is usually as follow:

1. Compute the setady-state thickness using the [thickness solver](#thickness_direct_solver)
2. Compute a cost function that measures the difference between the model and some observations, usually the observations are discrete along flow lines and the *Adjoint_CostDiscSolver* will be used.
3. Impose some constraint on the spacial variations of $H$. This can be done by penalising first spatial derivatives of $H$ with the *Adjoint_CostRegSolver* (This has to be done before the next step as it involves the variable $H$.)    
4. Compute the solution of the [adjoint linear system](#adjoint_linearsolver)
5. **Compute the gradient of your cost function with respect to your input parameters**

Note that the step 3 is usually crucial as it will impose spatial correlation in the retrieved ice thickness field and thus drives the results far from the observations.

### Keywords


#### Solver Section:

```
Solver *id*
  Equation = "DJDp"

  procedure = "ElmerIceSolvers" "AdjointThickness_GradientSolver"

 !# as for the direct solver the convection velocity 
 !#  can be provided here or in the body forces
  Flow Solution Name = String "...."

 !# name of the variable solution of the direct problem
  Thickness Solution Name = String ...

 !# name of the variable solution of the adjoint problem
  Adjoint Solution Name = String ...

 !# Compute the derivatives

 !# with respect to the velocity; results will be stored in the variable DJDUV
  ComputeDJDUV = Logical ... [default: False]
 !# Should we initialise the gradient to 0 here?
  Reset DJDUV = Logical ... [default: True]

 !# with respect to the top surface accumulation; 
 !#  results will be stored in the variable ComputeDJDsmbTop
  ComputeDJDsmbTop = Logical ... [default: False]
  Reset DJDsmbTop = Logical ... [default: True]
 
 !# with respect to the top surface accumulation;
 !#  results will be stored in the variable ComputeDJDsmbBot
  ComputeDJDsmbBot = Logical ... [default: False]
  Reset DJDsmbBot = Logical ... [default: True]

end


```

#### Body Forces:

```
!# Ms and Mb
Body Force *id*
  Top Surface Accumulation = Real ...
  Bottom Surface Accumulation = Real ...

  !# If <Flow Solution Name> is not provided in the solver section
  !#  the convection velocity can be read here:
  ! x-componenent
  Convection Velocity 1 = Real ...
  ! y-componenent if Convection Dimension = 2
  Convection Velocity 1 = Real ...
End
````

### Tests and Examples

- See examples for the [Mass Conservation methods](../../examples/Inverse_Methods)
