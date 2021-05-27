# Solver Shallow Ice Approximation (SIA)
## General Information
- **Solver Fortran File:** SIASolver.f90
- **Solver Name:** SIAVariable and SIASolver
- **Required Output Variable(s):** SIAFlow
- **Required Input Variable(s):** Depth, FreeSurfGrad1 and FreeSurfGrad2
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This SIA solver is not classical in the sense that the equations are not solved on a grid of dimension lower than the problem dimension itself. The geometry (H, B and S) is here given by the mesh. For a flow line problem, the mesh is a plane surface, and a volume for a 3D problem. Regarding this aspect, this solver is certainly not as efficient as a classical SIA solver. But, on the other hand, it works for an unstructured grid and non-constant viscosity. The SIA velocities and pressure can be used, for example, as initial conditions for the Stokes Solver. Contrary to the NS solver, the gravity must be orientated along the z-axis.

The SIA solver uses the same input parameters as the NS solver (Viscosity, Density, Viscosity Exponent, Flow BodyForce,…). Doesn't work with the built-in Glen's flow law (TODO). Need to run the [FlowDepth solver](./FlowDepth.md) first.

The basal velocities are given as Dirichlet BC on the bedrock surface. The [SSA Solver](./SSA.md) can be used for this purpose.

More information regarding how the SIA is solved using the SIASolver can be found [here](./siasolver.pdf).

## SIF contents
The required keywords in the SIF file for this solver are:

```
! Dummy solver just here to declare SIAFlow 
! as a true variable (not an exported variable)
! to allow access to previous values
Solver 2
  Equation = "SIA Variable"
  Procedure = File "ElmerIceSolvers" "SIAVariable"
  Variable = "SIAFlow"
  Variable DOFs = 4  ! 4 in 3D (u,v,w,p), 3 in 2D (u,v,p)
End

Solver 3
  Equation = "SIA"
  Procedure = File "ElmerIceSolvers" "SIASolver"
  Variable = -nooutput "SIAvar"
  Variable DOFs = 1

  Linear System Solver = "Direct"
  Linear System Direct Method = umfpack

  Steady State Convergence Tolerance = Real 1.0e-3
End

!!! bedrock
Boundary Condition 5
  Target Boundaries = 5
  SIAFlow 1 = Real 0.0e0
  SIAFlow 2 = Real 0.0e0
  SIAFlow 3 = Real 0.0e0
End

!!! free surface
Boundary Condition 6
  Target Boundaries = 6
  Save Line = Logical True
  SIAFlow 4 = real 0.0  !(p=0)
  Depth = real 0.0
End
```

## Examples
An example using the SIASolver applied to experiment A160 of ISMIP-HOM benchamrks can be found in [ELMER_TRUNK]/elmerice/Tests/SIA.
