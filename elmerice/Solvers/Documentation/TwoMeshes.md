# Solver TwoMeshes

## General Information

- **Solver Fortran File:** TwoMeshes.F90
- **Solver Name:** TwoMeshes
- **Required Output Variable(s):** Height
- **Required Input Variable(s):** List of Vars to interpolate (see example)
- **Optional Output Variable(s):** 
- **Optional Input Variable(s):** List of vars to nullify (set to zero)

## General Description
This is effectively the 2D version of CalvingRemesh, except it doesn't need to actually remesh, it simply deforms it instead.

## SIF contents

```
Solver 18

  Equation = MeshToMesh

  Procedure = "ElmerIceSolvers" "TwoMeshes"
  Variable = Height
  Exec Solver = "After Timestep"
!  Exec Solver = "Never"
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-8
  Linear System Preconditioning = ILU1
  Linear System ILUT Tolerance = 1.0e-3
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1
  Linear System Precondition Recompute = 1

!  Map Condition = Equals CalvingOccurs
!  Extension Amplitude = Real 0.5
  FS Bottom = Logical True
  FS Top = Logical True
  Nullify 1 = String "Long Mesh Update"
  Nullify 2 = String "Vert Mesh Update"
  Nullify 3 = String "Front Displacement"

  Variable 1 = Velocity 1
  Variable 2 = Velocity 2
  Variable 3 = Pressure
  Variable 4 = Stress 1
  Variable 5 = Stress 2
  Variable 6 = Stress 3
  Variable 7 = Stress 4
  Variable 8 = Depth
  Variable 9 = Height
  Variable 10 = StrainRate
  Variable 11 = Min Zs Bottom
  Variable 12 = Mesh Velocity 1
  Variable 13 = Mesh Velocity 2
  Variable 14 = weight
  Variable 15 = bedrock
  Variable 16 = fwater 1
  Variable 17 = fwater 2
  Variable 18 = Water Pressure
  Variable 19 = Distance

  Surface Variable 1 = String "Zs Top"
  Surface Variable 2 = String "Zs Bottom"
  Surface Variable 3 = String "GroundedMask"
End
```
