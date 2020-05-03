---
title: |
  | ElmerIce Documentation :
  | Adjoint inverse methods
author:
- F. Gillet-Chaulet
date: 
- \today
---

## Introduction 

This document provides documentation about solvers and methods implemented in Elmer/Ice
to solve optimisation problems with adjoint methods.

Currently these methods are restricted to diagnostic, i.e. steady state simulations, to optimise model parameters.

In general the workflow of an optimisation problem in Elmer/Ice is as follow:

1. Compute the solution of a problem with a solver. In general this will lead to solving the following linear system:  
$$ \boldsymbol{K(p)}.\boldsymbol{x} = \boldsymbol{F(p)} $$  
with $\boldsymbol{x}$ the solution depending on the input parameters $p$, $\boldsymbol{K}$ the stiffness matrix and
$\boldsymbol{F}$ the force vector. If the problem is non-linear this is not properly taken into account, however
if the solver is equiped with a Newton method, the gradients should be relatively accurate.
2. Compare the solution $\boldsymbol{x}$ with some observations, and evaluate a *cost function* $J(p)$.
3. Solve the adjoint of the linear system (step 1).
4. Compute the derivatives of $J$ with respect to $p$ from the solutions of the direct and adjoint problems in steps 1 and 3.
5. Add some regularisation term. If regularisation depends on $\boldsymbol{x}$ this should be done before step 3.
6. Search the set of parameters that minimises $J$.


[Chapter 1](#generic_solvers) provides documentation about generic solvers used for the validation, the optimisation (step 6)
and the solution of the adjoint of the linear system (step 3).

[Chapter 2](#cost_solvers) provides documentation for the solvers that can be used to compare the solution with observations (step 2). 
Some solvers are model specific.

[Chapter 3](#reg_solvers)  provides documentation for the solvers that can be used for regularisation (step 5). 
Some solvers are model specific.

The following chapters describe solvers specific to compute the gradients (4), and are model specific.

The adjoint codes have been derived by hand as usually the stucture of Elmer is too complex
to use  automatic differentiation softwares.

The solvers have been designed to be as generic as possible. A strength of Elmer is its modularity,
however it means that the configuration can be complex, and it is possible that when writting the code
whe have not taken into account all the possibilities offered by Elmer.

We have done the maximum to test the solvers and they should be accurate for standard simulations. However, it will
be to the user responsability to check that his configuarations leads to accurate gradients and smooth optimisation. 
A cost which is not decreasing is often the sign the the gradients are not accurates.

**Please contribute to improve this documentation and Elmer/Ice capabilities** by reporting errors or inaccuracies,
and contribute in developping new test cases and functionnalities.


