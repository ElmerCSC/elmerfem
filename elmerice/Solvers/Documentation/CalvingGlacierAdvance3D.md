# Solver/User Function GlacierAdvance3D
## General Information
- **Solver/User Function Fortran File:** CalvingGlacieradvance3D.F90
- **Solver/User Function Name:** Glacieradvance3D
- **Required Output Variable(s):** FrontAdvance
- **Required Input Variable(s):** 
- Normal Vector Variable Name = String "Front Normal Vector"
- Flow Solution Variable Name = String "Flow Solution"
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):**
- Melt Variable Name = String "FrontMelt"
- **Solver Keywords:**
- Ignore Velocity = Logical False. Permits frontal melt only
- Left Rail File Name = File ... location of left rail xy file
- Left Rail Number Nodes = Integer 12. The number of rail nodes in the above file.
- Right Rail File Name = File ... location of right rail xy file.
- Right Rail Number Nodes = Integer 12. The number of rail nodes in the above file.
- Rail Buffer = Real 0.1. The acceptable error that lateral nodes can deviate from the rails.
  
## General Description
This solver predicts the terminus advance based of the velocity field, a set of fjord rails and a melt field. The details of how this fits into the calving algorithm can be found [here](https://zenodo.org/records/10182710).

## Known Bugs and Limitations
Can cause minor artifical mass change. See calving algorithm document for more details.

## SIF Contents
The required keywords in the SIF file for this solver/USF are:

```
Solver n
  Equation = "Front Advance"
  Procedure =  File "ElmerIceSolvers" "GlacierAdvance3D"
  Exec Solver = "After Timestep"
  Save Exec When = String "After Timestep"
  Variable = -dofs 3 "FrontAdvance"

  Normal Vector Variable Name = String "Front Normal Vector"
  Flow Solution Variable Name = String "Flow Solution"
  Ignore Velocity = Logical False !Permits frontal melt only

  Left Rail Number Nodes = Integer 12
  Right Rail Number Nodes = Integer 12

  Rail Buffer = Real 0.1
End

Solver n+1 ! modify the mesh using front advance variable
  Equation = "Longitudinal Mesh Update"
  ! usually, the solver is executed only after the thermo-mechanical
  ! problem has obtained a solution on the time-level
  Procedure =  File "MeshSolve" "MeshSolver"
  Exec Solver = "After Timestep"
  Save Exec When = String "After Timestep"
  Solver Timing = Logical True

  Variable = Longitudinal Mesh Update
  Variable DOFs = 3
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations  = 500
  Linear System Preconditioning = ILU1
  Linear System Convergence Tolerance = 1.0e-12
  Linear System Abort Not Converged = False
  Nonlinear System Max Iterations = 1
  Nonlinear System Convergence Tolerance = 1.0e-06
  First Time Non-Zero = Logical True
End

Boundary Condition 1 ! front
  Longitudinal Mesh Update 1 = Equals FrontAdvance 1
  Longitudinal Mesh Update 2 = Equals FrontAdvance 2
  Longitudinal Mesh Update 3 = Equals FrontAdvance 3
End

Boundary Condition 2 ! sidewalls
  Longitudinal Mesh Update 1 = Equals FrontAdvance 1
  Longitudinal Mesh Update 2 = Equals FrontAdvance 2
  Longitudinal Mesh Update 3 = Equals FrontAdvance 3
End
```

To restrict the mesh update to only a small section of the domain add:
```
Body Force
  Longitudinal Mesh Update 1 = Real 0.0
  Longitudinal Mesh Update 2 = Real 0.0
  Longitudinal Mesh Update 3 = Real 0.0

  Longitudinal Mesh Update 1 Condition = Variable Distance
    Real MATC "tx-1500"
  Longitudinal Mesh Update 2 Condition = Variable Distance
    Real MATC "tx-1500"
  Longitudinal Mesh Update 3 Condition = Variable Distance
    Real MATC "tx-1500"
End
```
## Examples
An example using the full calving algorithm can be found here [ELMER_TRUNK]/elmerice/Tests/Calving3D_lset

## References
Iain Wheel PhD thesis - DOI: https://doi.org/10.17630/sta/611.
Full calving algorithm detailed [here](https://zenodo.org/records/10182710).
