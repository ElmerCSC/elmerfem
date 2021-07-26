## Optimisation with m1qn3

**Module name**: Optimize_m1qn3Parallel  
**Module subroutines**: Optimize_m1qn3Parallel  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)  
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 03.05.2020


### Introduction

[M1QN3](https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.html) 
is a routine designed to solve large scale unconstrained minimization problems. 
It implements a limited memory quasi-Newton algorithm (L-BFGS).

For details on the routine see the [documentation](https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.pdf).

This solver calls M1QN3 in *reverse* communication, *i.e.* the optimisation loop is controlled 
by running steady-state iterations.
Each iteration  computes the cost function $J$ and its gradient $g$ for the current state then pass these informations
to M1QN3 which provive a new estimate of the state.

### Specific comments

- Contrary to what his name suggests this solver now works for serial and parallel simulations.

- The *optimised variable* can be a vector in Elmer, i.e. a variable with DOFs>1. This property can be used
for example if we want to optimise both the friction coefficient and the viscosity with the SSA. In this case the
gradient should also be a vector with the same DOFs.

- The M1QN3 routine is *serial* so for parallel simulations, process 0 gather the state and gradients from all other partitions
and call M1QN3 in serial; this could become a bottle-neck for very large simulations.

- The keyword *M1QN3 df1* is the expected relative decrease of $J$ at the first iteration; it is used to scale the 
step-size at the first iteration. In general this number should not be too small; However if you start very close from
the minimum setting a too high value might require several simulations in the line search before $J$ can be effectively decreased.
We did not implement the possibility of a *warm restart*, so in gerenal we do not advise to restart an optimisation; however;
in this case, it might be beneficial to decrease this value to save some computing resssources in the first iterations.

- By default we use the  **Euclidean inner product**, which is also the inner product used to compute the gradients. 
In general the gradient are the results of some integration over the elements, so larger elements will lead to larger 
gradients. Experience has shown that this is a sensible choice with unstructured meshes where large elements correspond
to small velocities, so that we obtained a good compromise during the optimisation between areas of low velocities - low errors and 
area of high velocities - high errors. However, for  applications with unstructured meshes where the mesh size
is not related to the magnitude of the errors, it may be beneficial to change the inner product and scale the gradients by the mesh size.
This is automatically achieved when setting the keyword *Mesh Independent=TRUE*.
This will lead to smoother gradients, however for areas with small errors this will lead to small gradients and theses areas will not be affected by the optimisation. In this case it might be beneficial to have a cost function that also measures the relative error.  
*Note* that the change of inner product affects the *conditionning* of the minimisation, so in *perfect simulations* both solutions should lead to the same exact minimum....

- Outputs from M1QN3 are saved in the *M1QN3 OutputFile*. The  verbosity level is controlled with the keyword *impres*.

- The keyword *M1QN3 ndz* corresponds to the number of updates *m* used to approximate the Hessian. The memory requirement increases
with this parameter. From M1QN3 documentation, a good compromise should be obtained by taking *m* between 5 and 10.

- The decrease of *J* in M1QN3 is enforced by the Wolfe line-search. 
For M1QN3 there is a new *iteration* when $J$ has been sufficiently decreased during the line search. 
So for a given *iteration* there could be  several *simulations*. We don't make this disctinction in Elmer, so each
**elmer steady-state iteration** correspond to *1 simulation* for M1QN3.

- In general, if the gradient is accurate there should be a small number of *simulations* per *iteration*. If M1QN3
is blocked in a line search, this means that $J$ can not be sufficiently decreased in the current direction 
and this could be the sign of an in-accurate gradient. The minimal step-size for the line search is controlled
with the keywords *dxmin*.

- The **normal way** of stopping for M1QN3 is when the norm of the gradient has been decreased by a sufficient value
controlled with the *epsg* keyword.

- In general, fixing this keyword in advance might be difficult and it might be appropriate to give a small value and
fix the maximum number of simulations and iterations as the same number os maximum steady state iterations in the simulation section.

- The steady state iterations in Elmer are controlled with the keywords *Steady State Min Iterations* and *Steady State Max Iterations*
in the simulation section and with the keywords *"Steady State Convergence tolerance = Real ..."* in the different solvers. Setting this
last value to a too small value may stop the simulation too early. On the contrary if the criterion on *epsg* is reached by M1QN3 
the same state will be used for the next iterations and the *Steady State Convergence tolerance* should be satified stopping the Elmer Simulation.

### Keywords

In general this solver will be the last solver of an assimilation sequence. Below are the keywords:

```
Solver id
  Equation = "Optimize_m1qn3"
  procedure = "ElmerIceSolvers" "Optimize_m1qn3Parallel"

 !# name of the variable that contains value of *J* for the current state
  Cost Variable Name = String ""  ! [default: CostValue]
 !# name of the *state variable* x
  Optimized Variable Name = String "" ! [default: Beta]
 !# name of the gradient of *J* w.r.t. *x*:
  Gradient Variable Name = String ""  ! [default: DJDB]

 !# [Optionnal] : If provided,
 !   x will not be affected by the optimsation where mask.lt.0
   Optimisation Mask Variable = String ""

 !# [Optionnal] : compute and save the euclidean norm of *g*
  gradient Norm File = String ""

  !Note: it may be beneficial to set this to True, which scales
  !the gradient by 1/boundary_weights. With this set to false,
  !larger elements produce larger gradients.
  Mesh Independent = Logical False

! M1QN3 Parameters
!  see M1QN3 documentation for full explanation
  M1QN3 dxmin = Real ... ! [default: 1.e-10]
  M1QN3 epsg = Real  .... ! [default: 1.e-6]
  M1QN3 niter = Integer ... ! [default: 200]
  M1QN3 nsim = Integer ... ![default: 200]
  M1QN3 impres = Integer ... ![default: 5]
  M1QN3 DIS Mode = Logical ... ![default: False]
  M1QN3 df1 = Real ...   ![default: 0.2]
  M1QN3 normtype = String ... ![default: dfn]
  M1QN3 OutputFile = File "" ![default: M1QN3.out]
  M1QN3 ndz = Integer ... ![default: 5]
end
```

### Tests and Examples

