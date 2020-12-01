#Solving the Mass Conservation of Snow/Firn

##General Description
This page explains how to use the general AdvectionReactionSolver from the Elmer distribution to get the density evolution in case of a compressible material (snow/firn) under a given velocity field computed from the [Porous Solver](./PorousSolve.md). The AdvectionReactionSolver solves the general equation

*{{\partial A}/{\partial t}} + div (A u) + gamma A=sigma*

where *u* is the velocity vector. In the particular case of the mass conservation equation, one has therefore *gamma = 0* and *sigma = 0*. Solving for the true density (kg/m^3) or the relative density is equivalent (but limit values and Dirichlet boundary conditions have to be set accordingly).

Note 1: equation (4.1) in the Elmer Model Manual for the AdvectionReaction sover is not correct. The previous equation is the one implemented.

Note 2: Have a look to this [post](http://elmerfem.org/forum/viewtopic.php?f=7&t=3066&p=9570#p9570) on the Elmer Forum regarding the initialisation of both the DG primary and exported variables of the AdvectionReaction solver (see the example at the end of this page).

##SIF contents
The Solver section looks like:

```
Solver 8
  Equation = "AdvReact"
  Exec Solver = "After Timestep"
  Procedure = File "AdvectionReaction" "AdvectionReactionSolver"
  ! this is the DG variable, which is not part of the output
  Variable =  -nooutput "DGdens"
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
  Exported Variable 1 = Relative Density
  Exported Variable 1 DOFS = 1
End
```
The source in case of the mass conservation equation is 0

```
Body Force 1
 ...
  DGDens Source = Real 0.0 
 End
```
Initial and boundary conditions are then being set for the DG variable and not the exported one!

```
Initial Condition 1
  ...
  DGDens = Real 0.4
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
  ! relative density on the upper surface
  DGDens = Real 0.4
End
```
The Material section contains a zero reaction rate as well as the upper/lower limit for the DG variable

```
Material 1
 ..
 ! Relative density must stay < 1
 DGDens Upper Limit = Real 1.0

 ! a minimum relative density is recommended for the Porous solver 
 DGDens Lower Limit = Real 0.3

 !Reaction rate is equal to zero
 DGDens Gamma = Real 0.0
End
```
##Examples
A 1D example build from an analytical solution can be found in [ELMER_TRUNK]/elmerice/Tests/Density. In that case, the velocity and density are inversely proportional (u(z) = K/D(z)).

A 3D example coupling the [Porous Solver](./PorousSolve.md) and the calculation of the density field can be found in [ELMER_TRUNK]/elmerice/Tests/DGsolver.
