# Solver OutputStrainHeating
## General Information
- **Solver Fortran File:** OutputStrainHeating.f90
- **Solver Name:** getStrainHeating

## General Description
Contains subroutine for writing output for strain heating (deformational heat) when strain heating is taken into account by setting Friction Heat = Logical True in the Body Force Section. This heat source is taken into account or not independently from running this Solver.
Do not use this subroutine for calculation (not precise enough) of strain heating, it's only meant for output purposes to see Friction Heat in the results files.
In case of convergence problems, it is recommended to use the [DeformationalHeat](./DeformationalHeat.md) Solver instead of the keyword Friction Heat = Logical True.

## SIF contents
```
Solver 1
  Equation = "StrainHeating"
  Variable = String "StrainHeat"
  Variable DOFs = 1
  Procedure = File "ElmerIceSolvers" "getStrainHeating"
  Nonlinear System Max Iterations = 1
End

Body Force 1
....
  Friction Heat = Logical True
....
End
```
