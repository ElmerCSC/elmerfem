## Adjoint Stokes : Slip Coefficient Gradient Solver {#stokes_gradientBeta_solver}

**Module name**: AdjointStokes_GradientBetaSolver.F90  
**Module subroutines**: AdjointStokes_GradientBetaSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)    
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 17/06/2020  

**Required input variables:**
   
 - Variable solution of the direct problem (Flow Solution)
 - The adjoint variable computed from the [adjoint linear system](#adjoint_linearsolver)

**output variables:**

 - Nodal derivatives with respect to the Slip Coefficient

### Introduction

This is the adjoint code corresponding to the application of the slip coefficients in the Stokes Solver.


The sequence in the .sif is usually as follow:

1. Compute the velocity using the Stokes solver
2. Compute a cost function that measure the diffrence between the model and some observations
3. Compute the solution of the [adjoint linear system](#adjoint_linearsolver)
4. **Compute the gradient of your cost function with respect to your slip coefficient**


If a change of variable is used for the input slip coefficient, the derivative of
the friction coefficient with respect to your input variable can directly be provided
in the *Boundary Condition* section. 

Be carefull:  
  - **This solver must be executed on the bottom boundary where the slip condition is applied**   
  - **The NormalTangential system has to be used**  
  - **In 3D Slip Coefficient 2 and 3 have to be equal**  
  - **It will initialise the derivative to 0 by default.** Use *Reset Gradient Variable=Logical False* if the gradient should be added to the values computed in another solver.   

### Keywords

#### Solver section

```
  Solver *id*
   
  procedure = "ElmerIceSolvers" "AdjointStokes_GradientBetaSolver"


  !# name of the variable solution of the direct problem 
   Flow Solution Name = String ... *[default:ssavelocity]*

  !# name of the variable solution of the adjoint problem 
   Adjoint Solution Name = String ... *[default:adjoint]*

  !# name of the variable that contains the derivative with respect to the slip coefficient
   Gradient Variable Name = String .... *[default: DJDB]*

  !# Should the gradient be initialised to 0
  Reset Gradient Variable = Logical ... *[default: True]*

  End

```
#### Boundary Conditionds

By default the output for the gradient is the derivative with respect to the slip coefficient.
If a change of variable is used, you can directly provide the derivative here:

e.g.
```
Boundary Condition *id*
 ! This boundary must be a body to execute the solver
  Body Id = *BdId*

  Normal-tangential Velocity = True
  Velocity 1 = Real 0.0e0

  Slip Coefficient 2 = Variable alpha
    REAL procedure "ElmerIceUSF" "TenPowerA"

  Slip Coefficient 3 = Variable alpha
    REAL procedure "ElmerIceUSF" "TenPowerA"

  Slip Coefficient derivative = Variable alpha
    REAL procedure "ElmerIceUSF" "TenPowerA_d"
End

```

### Tests and Examples

- See examples for the [Stokes inverse methods](../../examples/Inverse_Methods)
