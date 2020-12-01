#ForceToStress Solver
##General Information
- **Solver Fortran File:** ForceToStress.f90
- **Solver Name:** ForceToStress
- **Required Output Variable(s):** Stress (user defined)
- **Required Input Variable(s):** Force
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

##General Description
For a given boundary, this solver computes the nodal stress equivalent to a given nodal force. Force here means the result of a variable, named stress, which is of the same units than the force but per unit-area in 3D or unit-length in 2D. This can be used also to infer the flux from the value of a debit or to infer a tangential stress from the value of a tangential force.

##SIF contents
In the SIF example below, the normal stress on a boundary is inferred from the 3rd component of the Stokes residual.

```
Solver 1
  Equation = "Navier-Stokes"
  Stabilization Method = String Stabilized
  Flow Model = Stokes
  ...

  Exported Variable 1 = Flow Solution Loads[Fx:1 Fy:1 Fz:1 CEQ Residual:1 ]
  Calculate Loads = Logical True
End

Solver 2
    Equation = "ForceToStress"
    Procedure = File "ElmerIceSolvers" "ForceToStress"
    Variable = String "Stress"
    Variable DOFs = 1
    
    Force Variable Name = String "Fz"

    Linear System Solver = Direct
    Linear System Direct Method = umfpack
 End

! Solve this for body Id 2 (=boundary 3 here)
Equation 2
  Active Solvers(1) = 2
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Boundary Condition 3
  Target Boundaries = 5
  Body Id = 2
  ...
End
```

##Examples
In the example found in [ELMER_TRUNK]/elmerice/Tests/ForceToStress, a pressure applied on a boundary is first integrated to get nodal force using the [GetHydrostaticLoad Solver](./GetHydrostaticLoads.md), and then the pressure is recovered using the ForceToStress Solver.

2 more tests were added on 15-12-2018:

[ELMER_TRUNK]/elmerice/Tests/ForceToStress_periodic: Test periodic conditions

[ELMER_TRUNK]/elmerice/Tests/ForceToStress_parallel: Test parallel solver
