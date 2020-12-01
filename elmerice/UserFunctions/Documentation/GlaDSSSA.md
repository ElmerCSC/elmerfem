# User Function UFS_CouplingGlaDS_SSA
## General Information
- **USF Fortran File:** USF_CouplingGlaDS_SSA.F90
- **USF Name:** HorizontalVelo and OverburdenPressure
- **Required Input Variable(s):** SSAVelocity, Zb and Zs
- **Optional Input Variable(s):** None

## General Description
The aim of this user function is to pass information from the SSA solver to the GlaDS solver.

## SIF contents
The required keywords in the SIF file for these user functions are:

```
Material 1
  ...
  Sliding Velocity = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "HorizontalVelo"
  
  Ice Normal Stress = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "OverburdenPressure"
  ...
End
```

## Examples
An example can be found in [ELMER_TRUNK]/elmerice/Tests/GlaDS_SSA.
