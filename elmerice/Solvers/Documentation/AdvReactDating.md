# Solving the age equation
## General Description
This page explains how to use the general AdvectionReactionSolver from the Elmer distribution to get an Age/Depth relation:

*{{\partial A}/{\partial t}} + u . grad(A) =1*

The AdvectionReactionSolver solves the general equation

*{{\partial A}/{\partial t}} + div (A u) + gamma A=sigma*

In the particular case of the age equation, one has *gamma = -epsilon_m* and *sigma = 1*. Obviously *gamma = 0* for an incompressible flow. For a compressible flow (snow/firn rheology), the trace of the strain-rate tensor *epsilon_m* should be evaluated using the [ComputeStrainRate Solver](./ComputeStrainRate.md) or directly from the [Porous Solver](./PorousSolve.md).

## SIF contents
The Solver section looks like:

```
Solver 8
  Equation = "AdvReact"
  Exec Solver = "After Timestep"
  Procedure = File "AdvectionReaction" "AdvectionReactionSolver"
  ! this is the DG variable, which is not part of the output
  Variable =  -nooutput "DGAge"
  ! this tells that the solver is run on DG mesh
  Discontinuous Galerkin = Logical True
  ! the solver can account for upper and lower limits of the variable
  ! imposed by formulation of an variational inequality (VI)
  ! next line switches the VI to be accounted for
  Limit Solution = Logical True

  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations  = 1000
  Linear System Preconditioning = ILU1
  Linear System Convergence Tolerance = 1.0e-06
  ! Variational inequality makes it a non-linear problem
  Nonlinear System Max Iterations = 40
  Nonlinear System Min Iterations = 2
  Nonlinear System Convergence Tolerance = 1.0e-04

  ! This is the variable that is used to interpolate
  ! the DG solution to the regular FEM mesh in order
  ! to get a correct output
  Exported Variable 1 = Age
  Exported Variable 1 DOFS = 1
End
```
The source in case of the age/depth equation is 1 (da/dt + u da/dx + … = 1)

```
Body Force 1
 ...
  DGAge Source = Real 1.0  ! result in years
 End
```
Initial and boundary conditions are then being set for the DG variable and not the exported one!

```
Initial Condition 1
  ...
  DGAge = Real 0.0
End

! only Dirichlet BC can be set
! the solver automatically uses this
! condition only on inflow boundaries
! outflow boundaries are ignored
Boundary Condition 2
  Name = "surf"
  Target Boundaries = 2
  Body ID = 2
  ...
  ! fresh snow or ice
  DGAge = Real 0.0
End

! a suggestion for a trick at no-slip BC's to
! keep the age in bounds
!----------------------------------------------
Boundary Condition 1
   Name = "bedrock"
   Normal-Tangential Velocity = True
  ! normal outflow velocity set to a very small value
  ! in order to avoid infinitely large age values
  Velocity 1 = Real 1.0e-06 (a micro-meter per year)
  ! no-slip tangential
  !-------------------
  Velocity 2 = Real 0.0
  Velocity 3 = Real 0.0
  .
  .
  .
End
```
Just make sure that a free surface can have both, inflow and outflow. NB: If you have a model for re-freezing at the bedrock, you also should set a DGAge = Real 0.0 on this boundary.

The Material section contains the reaction rate as well as the upper/lower limit for the DG variable

```
Material 1
 ..
 ! this has to be sometimes set in diagnostic
 ! runs to avoid infinite [actually O(1/machine precision)]
 ! ages at stagnation points in the flow field
 DGAge Upper Limit = Real 100000

 ! well, no ice from the future
 ! anyway, if you have a well behaving velocity field
 ! you should not be in need of that constraint
 DGAge Lower Limit = Real 0

 !this would be a reaction rate,
 ! is equal to -tr(Eij) (minus trace of the strain-rate) 
 ! in the case of an incompressible material (ice), it is then 0
 DGAge Gamma = Real 0.0
End
```
## Example
An example using the AdvectionReaction solver to solve the age equation can be found in [ELMER_TRUNK]/elmerice/examples/Test_Dating. For a steady 1D flow, the dating equation simplifies to dA/dx = 1/u and solutions are of the form A = ln(Kv) with u = v/v'.
Two tests can be also found in [ELMER_TRUNK]/elmerice/Tests/Dating and [ELMER_TRUNK]/elmerice/Tests/Density.
