## AdjointSSA : Flux divergence cost solver {#ssa_fluxdiv_cost}

**Module name**: AdjointSSA_CostFluxDivSolver.F90  
**Module subroutines**: AdjointSSA_CostFluxDivSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)    
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 28/04/2020  

**Required input variables:**
   
 - Variable solution of the direct problem ([SSA solver](#ssa_direct_solver)) (SSAVelocity *[hard coded]*)
 - top surface elevation : Zs *[hard coded]*
 - bottom surface elevation : Zb *[hard coded]*

**output variables: **

- CostValue [computed cost can be a global variable]
- Velocityb [derivative with respect to the ssavelocity]
- *[Optional]* DJDZb [derivative with respect to Zb]
- *[Optional]* DJDZs [derivative with respect to Zs]


### Introduction

Compute a cost function that compares the *diagnostic* model thickness rate of change
$$ dhdt^{mod}=\dfrac{\textrm{d}H}{\textrm{d}t} = -\textrm{div}(uH) + M_t + M_b$$
with the observed $dhdt^{obs}$ as
$$J_{fdiv}= \int_\Omega 0.5 (dhdt^{obs} - dhdt^{mod} )^2 d\Omega$$

  - $\textrm{div}(uH)$ is the diagnostic flux divergence
  - $M_t$ the top surface mass balance (>0 for accumulation)
  - $M_b$ the bottom surface mass balance (>0 for accumulation)

In general flux divergences anomalies are due to uncertainties in the bedrock topography and others model assumptions,
so that adjusting only the friction parameter will not allow to perfectly fit the observations.   
But adding this cost function in addition to the direct comparison with observed surface velocities may allow to
reduce initial flux divergence anomalies. Putting too much weight to this cost function may tend to decrease the
magnitude of the velocity at the expense of observations.

As this cost function depends on the solution of the SSA, it has to be run before the adjoint of the linear system.

The sequence in the .sif will usually be as follow:

1. Compute the velocity using the [SSA solver](#ssa_direct_solver)
2. Compute a cost function that measure the diffrence between the model velocities and some observation
3. **Compute the Flux divergence cost**
4. Compute the solution of the [adjoint linear system](#adjoint_linearsolver)
5. Compute the gradient of your cost functions with respect to your input parameters ([SSA Gradient Solver](#ssa_gradient_solver))


### Keywords

#### Solver section

```
  Solver *id*
   Equation = "Cost_DHDT"
   
   procedure = "ElmerIceSolvers" "AdjointSSA_CostFluxDivSolver"

  !# Should we initialise Cost and velocityb to 0 here
  !# Set to false is a cost function has been computed in a previous solver
   Reset Cost Value = Logical [default: true]

  !# Name of Cost Variable
   Cost Variable Name = String "CostValue"  

  !# Multiplicative scaling parameter for the cost
   Lambda= Real ....

  !# Should we compute gradient with respect to bottom surface elevation
  !# hard coded initialization to zero in this solver
   Compute DJDZb = Logical .... [default : True]
  !# Should we compute gradient with respect to top surface elevation
  !# hard coded initialization to zero in this solver
   Compute DJDZs = Logical .... [default : True]

 !# save the cost as a function of iterations (iterations,Cost,rms=sqrt(2*Cost/area))
   Cost Filename = File ...

  End

```

#### Body Force section

Provide values for the mass balance and observation in the body force section.

```
  BodyForce *id*
 
   Top Surface Accumulation = Real ....
   Bottom Surface Accumulation = Real ....

   Observed dhdt = Real ...
  End
```

You can use keywords related to passive elements if you want to skip the evaluation of the 
cost function for some elements. By default the name of the solver varaible is *[Equation_name]*_var.

```
  BodyForce *id*

  Cost_DHDT_var passive = Real ...

End
```

### Tests and Examples

