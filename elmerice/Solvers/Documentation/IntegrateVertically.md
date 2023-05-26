# IntegrateVertically Solver
## General Information
- **Solver Fortran File:** IntegrateVertically.f90
- **Solver Name:** IntegrateVertically
- **Required Output Variable(s):** given by the variable name
- **Required Input Variable(s):** Name of the variable to be integrated in Integrated Variable Name, Depth or Height
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This solver computes the depth integrated value of a variable (sum_zb^zs D dz) or the mean value (1/H sum_zb^zs D dz). The output of resulting integrated value or mean can be done at the upper surface or the lower surface. If the mean value is computed (Compute Mean = Logical True), the depth (if On Surface = Logical False) or the height (if On Surface = Logical True) has to be calculated, too. If the integrated variable is calculated on the upper surface (On Surface = Logical True), a Dirichlet BC has to be given at the bottom surface, and vice versa.

## SIF contents
```
Solver 2
  Equation = "IntegrateVertically"

   Procedure = File "ElmerIceSolvers" "IntegrateVertically"
   Variable = String "Mean Var"
   Variable DOFs = 1

   Exported Variable 1 = String "VarToBeIntegrated"
   Exported Variable 1 DOFs = 1

   ! We want it computed on the bed
   On Surface = Logical False
   ! We want the mean value
   ! We then need the Depth
   Compute Mean = Logical True
   Integrated Variable Name = String "VarToBeIntegrated"

   Linear System Solver = "Direct"
   Linear System Direct Method = umfpack
End

!!! free surface
Boundary Condition 2
  Target Boundaries = 6
  Mean Var = Real 0.0
  Depth = Real 0.0
End
```

## Examples
An example using the IntegrateVertically Solver can be found in [ELMER_TRUNK]/elmerice/Tests/IntegrateVertically.
