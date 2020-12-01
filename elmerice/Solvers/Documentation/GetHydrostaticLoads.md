#GetHydrostaticLoads Solver
##General Information
- **Solver Fortran File:** GetHydrostaticLoads.f90
- **Solver Name:** GetHydrostaticLoads
- **Required Output Variable(s):** Fw (user defined DOFs=2 in 2D, 3 in 3D)
- **Required Input Variable(s):** External Pressure
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

##General Description
For a given boundary, this solver computes the equivalent nodal force induced by a pressure on the surface (read from the keyword External Pressure). The result is a force vector (Fx,Fy) in 2D and (Fx,Fy,Fz) in 3D.

##Known bugs
2018-12-15: Periodic boundary conditions fixed

2013-03-07: GetHydrostaticloads has an issue for periodic boundary conditions: the force is half what it should be on the periodic boundary whereas the residual account for the fact that the domain is periodic. TODO: fix this!

##SIF contents
```
Solver 1
  Equation = Fw
  Variable = Fw[Fx:1  Fy:1 Fz:1]
  Variable DOFs = 3
  Procedure = "ElmerIceSolvers" "GetHydrostaticLoads"
End

! Solve this for body Id 2 (=boundary 3 here)
Equation 2
  Active Solvers(1) = 1
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Boundary Condition 3
  Target Boundaries = 5
  Body Id = 2
  External Pressure = Real 0.1
End
```

##Examples
In the example found in [ELMER_TRUNK]/elmerice/Tests/ForceToStress, a pressure applied on a boundary is first integrated to get nodal force using the GetHydrostaticLoads solver, and then the pressure is recovered using the [ForceToStress Solver](./ForceToStress.md).
