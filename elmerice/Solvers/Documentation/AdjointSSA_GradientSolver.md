## AdjointSSA : Gradient Solver {#ssa_gradient_solver}

**Module name**: AdjointSSA_GradientSolver.F90  
**Module subroutines**: AdjointSSA_GradientSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)    
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 28/04/2020  

**Required input variables:**
   
 - Variable solution of the direct problem ([SSA solver](#ssa_direct_solver))
 - The adjoint variable computed from the [adjoint linear system](#adjoint_linearsolver)

**output variables: [Optional]**

 - Nodal derivatives with respect to: 

    - SSA Friction Parameter [set in Material]
    - SSA mean viscosity  [set in Material]
    - SSA mean density    [set in Material]
    - nodal top surface elevation [Zs variable]
    - nodal bottom surface elevation [Zb variable]


### Introduction

This is the adjoint code of the [SSA solver](#ssa_direct_solver).

Computes the derivatives SSA with respect to the input parameters.

The sequence in the .sif is usually as follow:

1. Compute the velocity using the [SSA solver](#ssa_direct_solver)
2. Compute a cost function that measure the diffrence between the model and some observations
3. Compute the solution of the [adjoint linear system](#adjoint_linearsolver)
4. **Compute the gradient of your cost function with respect to your input parameters**


If a change of variable is used for the input friction coefficient or viscosity, the derivative of
the friction coefficient or viscosity with respect to your input variable can directly be provided
in the Material section. 


### Keywords

#### Solver section

```
  Solver *id*
   
  procedure = "ElmerIceSolvers" "AdjointSSA_GradientSolver"

  !# name of the variable solution of the direct problem 
   Flow Solution Name = String ... *[default:ssavelocity]*

  !# name of the variable solution of the adjoint problem 
   Adjoint Solution Name = String ... *[default:adjoint]*

  !# Should we compute the gradient w.r.t. **SSA Friction Parameter** [Material]
   Compute DJDBeta = Logical ....
  !# Name of the variable that will contain the derivative
   DJDBeta Name = String .... *[default:DJDBeta]*
  !# Should we initialise the gradient to 0 here?
   Reset DJDBeta = Logical  ....

  !# Same for **SSA Mean Viscosity** [Material]
  !# Var Name hard coded DJDEta
   Compute DJDEta = Logical ...
   Reset DJDEta = Logical ...

  !# Same for **SSA Mean Density** [Material]
  !# Var Name hard coded DJDRho
   Compute DJDRho = Logical ...
   Reset DJDRho = Logical ...

  !# Same for nodal Zs [Zs Variable]
  !# Var Name hard coded DJDZs
   Compute DJDZs = Logical ...
   Reset DJDZs = Logical ...

  !# Same for nodal Zb [Zb Variable]
  !# Var Name hard coded DJDZb
   Compute DJDZb = Logical ...
   Reset DJDZb = Logical ...

  End

```
#### Material Properties:

By default the output for DJDBeta and DJDEta is the derivative with respect to the values 
prescribed in the material properties. If a change of variable is used, you can directly provide the derivative here:

e.g.
```
 Material *id*

  SSA Friction Parameter = Variable alpha
    REAL MATC "10^tx"

  SSA Friction Parameter Derivative = Variable alpha
    REAL MATC "ln(10)*10^tx"

  SSA Mean Viscosity = Variable alpha
    REAL MATC "tx*tx"

  SSA Mean Viscosity Derivative = Variable alpha
    REAL MATC "2*tx"

 End
```
Or using ElmerIce build-in user functions:

```
 Material *id*

 SSA Friction Parameter = Variable alpha
   REAL procedure "ElmerIceUSF" "TenPowerA"

 SSA Friction Parameter Derivative = Variable alpha
   REAL procedure "ElmerIceUSF" "TenPowerA_d"

 
 SSA Mean Viscosity = Variable alpha
    REAL procedure "ElmerIceUSF" "Asquare"

  SSA Mean Viscosity Derivative = Variable alpha
    REAL procedure "ElmerIceUSF" "Asquare_d"

 End

```

### Tests and Examples

- See examples for the [SSA inverse methods](../../examples/Inverse_Methods)
