# IntegratedVelocity Solver
## General Information
- **Solver Fortran File:** IntegratedVelocity.f90
- **Solver Name:** IntegratedVelocity
- **Required Output Variable(s):** Integrated Velocity (DIM-1)
- **Required Input Variable(s):** A Flow Solution (given in Flow Solution Name), Depth or Height
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This solver compute the mean vertically integrated horizontal velocities, either on the bedrock surface or on the upper surface. If the keyword Compute Flux is set to True, the horizontal fluxes are computed else the mean horizontal velocities are computed (default).

## SIF contents
```
Solver 3
  Equation = "Integrated Velocity"
  Procedure = File "ElmerIceSolvers" "IntegratedVelocity"
  Variable = -nooutput "varVelo"
  Variable DOFs = 1

  Exported Variable 1 = String "Integrated Velocity"
  Exported Variable 1 DOFs = 2  ! 1 in 2D, 2 in 3D

  Flow Solver Name = String SSAFlow
  Compute Flux = Logical False
  On Surface = Logical True

  Linear System Solver = Direct
  Linear System Direct Method = umfpack
End

!Bedrock BC
Boundary Condition 1
  Target Boundaries = 1

  Integrated Velocity 1 = Equals SSAFlow 1
  Integrated Velocity 2 = Equals SSAFlow 2
  Height = Real 0.0
End
```

## Examples
An example using the IntegratedVelocity Solver can be found in [ELMER_TRUNK]/elmerice/Tests/IntegratedVelocity.
