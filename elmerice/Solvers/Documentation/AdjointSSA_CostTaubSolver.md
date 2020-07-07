## AdjointSSA : Basal drag regularisation solver {#ssa_Taub_cost}

**Module name**: AdjointSSA_CostTaubSolver.F90  
**Module subroutines**: AdjointSSA_CostTaubSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)    
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 28/04/2020  

**Required input variables:**
   
 - Variable solution of the direct problem ([SSA solver](#ssa_direct_solver)) (SSAVelocity *[hard coded]*)

**output variables:**

- Velocityb [*hard coded* derivative with respect to the ssavelocity]
- DJDBEta   [derivative with respect to the frcition parameter]
- CostValue [computed cost (can be a global variable)]


### Introduction

Instead of penalising the first spatial derivatives of the friction parameter using the [Adjoint_CostReg Solver](#adjoint_CostReg), this solver penalises the spatial derivatives of the basal drag $\tau_b$.  
$\tau_b$ is computed at the nodes from the parameters prescribed for the [SSA solver](#ssa_direct_solver).  

The cost function then writes:  
$$J_{\tau_b}=\int_\Omega 0.5 (\dfrac{\textrm{d}\tau_b}{\textrm{d} x})^2  \textrm{d} \Omega $$

If the friction is set to 0 automatically in the [SSA solver](#ssa_direct_solver) using the keyword 
*Sub-Element GL parameterization = logical True*, $\tau_b$ as computed here will be erroneous. 
If there is **floating elements**, you should skip the evaluation of this cost function in the floating elements
using the keywords related to passive elements, as, by definition, $\tau_b=0$ in floating elements.

As this cost function depends on the solution of the SSA, it has to be run before the adjoint of the linear system.

The sequence in the .sif will usually be as follow:

1. Compute the velocity using the [SSA solver](#ssa_direct_solver)
2. Compute a cost function that measure the diffrence between the model velocities and some observation
3. **Compute the regularisation term for $\tau_b$**
4. Compute the solution of the [adjoint linear system](#adjoint_linearsolver)
5. Compute the gradient of your cost functions with respect to your input parameters ([SSA Gradient Solver](#ssa_gradient_solver))

### Limitations and possible improvments

- Evaluation of $\tau_b$ do not check if the element is floating; see above.

- Similarly to the regularisation solver, a cost function that measures the difference between $\tau_b$
an a *prior* estimate could be implemented.


### Keywords

#### Solver section

```
  Solver *id*

    Equation = "DJDTau_Reg"
    procedure = "ElmerIceSolvers" "AdjointSSA_CostTaubSolver"

   !# Should we initialise cost to 0 here
    ! in general will be FALSE if a cost function has been 
    ! computed in a previous solver
     Reset Cost Value = Logical ...
  
   !# name of the cost file for output
     Cost Filename = File "CostReg_$name$.dat"

   !# name of the cost variable (can be a global variable)
     Cost Variable Name = String "CostValue"

   !# name of the variable that contains derivative with respect
   !  to the friction parameter
     DJDBeta Name = String  ....
   !# Initialise DJDBeta to 0 here? Has to be TRUE if this 
   !   solver is the first in the sequence to compute this variable
     Reset DJDBeta = Logical ....

   !# Multiplicative scaling parameter for the cost function
     Lambda = Real ....

  End

```

#### Material properties

By default the output for DJDBeta is the derivative with respect to the *SSA Friction Parameter*
prescribed in the material properties. 
If a change of variable is used, you can directly provide the derivative here:

e.g.
```
 Material *id*

  SSA Friction Parameter = Variable alpha
    REAL MATC "10^tx"

  SSA Friction Parameter Derivative = Variable alpha
    REAL MATC "ln(10)*10^tx"

 End
```
Or using ElmerIce build-in user functions:

```
 Material *id*

 SSA Friction Parameter = Variable alpha
   REAL procedure "ElmerIceUSF" "TenPowerA"

 SSA Friction Parameter Derivative = Variable alpha
   REAL procedure "ElmerIceUSF" "TenPowerA_d"

 End

```

#### Body Force section

You can use keywords related to passive elements if you want to skip the evaluation of the 
cost function for some elements, for example in floating elements.

By default the name of the solver varaible is *[Equation_name]*_var.

```
  BodyForce *id*

   DJDTau_Reg_var passive = Real ...

  End
```

### Tests and Examples

- See examples for the [SSA inverse methods](../../examples/Inverse_Methods)
