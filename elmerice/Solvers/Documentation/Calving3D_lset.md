# Solver/User Function Calving3D_lset
## General Information
- **Solver/User Function Fortran File:** Calving3D_lset.F90
- **Solver/User Function Name:** Calving3D_lset
- **Required Output Variable(s):** Calving Lset
- **Required Input Variable(s):** Stress for the auxiliary solver ProjectCalving.
  - Calving Stress Variable Name = String "Stress"
- **Optional Output Variable(s):** 
- **Optional Input Variable(s):** Binary crevasse field in he auxiliary solver ProjectCalving
  - Save Crevasse Field = Logical True
- **Solver Keywords:** 
Keywords in the solver:
- Calving Search Distance = Real 3000.0. The distance upstream crevasses are searched for.
- PlaneMesh Grid Size = Real 40.0. The element size of the 2D crevasse map.
- Lateral Calving Margins = Logical True. Account for the lateral margins when validating crevasses. See calving algorithm overview for more details.
- Crevasse Penetration = Real 1.0. Define the proportion of the ice column that needs to be crevassed to induce calving (value between 0-1).
- Project Calving Equation Name = String "CalvingProjection". Point to the auxiliary solver.
- Isosurface Equation Name = String "Calving Isosurface". Point to the auxiliary solver.
- Front Advance Solver = String "Front Advance". Point to the FrontAdvance solver. This is required so the fjord margins can be accounted for.
  
## General Description
This solver produces a level set or signed distance variable where the zero contour indicates the new terminus using the crevasse depth calving law. An overview of how the calving 3D solvers fit together is described in Calving.md. Further details of the new calving algorithm can be found [here](https://zenodo.org/records/10182710).

### Additional subsection 1 {e.g. Geographical Restriction}
Requires two auxiliary solvers as outlined in the example below. ProjectCalving and Calving Isosurface. Additionally, to save the 2D crevasse mesh a extra call of the result output solver is required.

## Known Bugs and Limitations
None known so far.

## SIF Contents
The required keywords in the SIF file for this solver/USF are:

```
Solver n
  Equation = "3D Calving Lset"
  Exec Solver = "After Timestep"
  Procedure = "ElmerIceSolvers" "Find_Calving3D_Lset"
  Solver Timing = Logical True

  Variable = String "Calving Lset"
  Variable DOFs = 1

  Exported Variable 1 = -dofs 1 "CIndex"
  Exported Variable 1 DOFs = 1

  Calving Search Distance = Real 3000.0
  PlaneMesh Grid Size = Real 40.0

  Lateral Calving Margins = Logical True

  Crevasse Penetration = Real 1.00

  Project Calving Equation Name = String "CalvingProjection"
  Isosurface Equation Name = String "Calving Isosurface"
  Front Advance Solver = String "Front Advance"
End

Solver n+1
  Equation = "CalvingResultOutput"
  Procedure = File "ResultOutputSolve2" "ResultOutputSolver"
  Exec Solver = "Never" !auxiliary solver called by Calving3D

  Output File Name  = "plane_$namerun"_" !"
  Vtu Format = logical true
  Binary Output = True
  Single Precision = True
  Save Geometry IDs = True

!  Scalar Field 1 = "ave_cindex"
End

Solver n+2
  Equation = "CalvingProjection"
  Procedure = File "ElmerIceSolvers" "ProjectCalving"
  Exec Solver = "Never" !auxiliary solver called by Calving3D
  Solver Timing = Logical True

  Basal Crevasse Model = Logical True
  Surface Crevasse Model = Logical True

  Calving Stress Variable Name = String "Stress"
  Plane Permutation(3) = Integer 1 3 2
  Volume Permutation(3) = Integer 1 3 2
End

Solver n+3
  Exec Solver = "Never" !auxiliary solver called by Calving3D
  Equation = "Calving Isosurface"
  Procedure = File "Isosurface" "IsosurfaceSolver"

  Isosurface Variable = String "ave_cindex"
  Isosurface Value = Real 0.0
  Isosurface Interpolant 1 = String "isoline id"
End
```
## Examples
An example using the full calving algorithm can be found here [ELMER_TRUNK]/elmerice/Tests/Calving3D_lset

## References
Iain Wheel PhD thesis - DOI: https://doi.org/10.17630/sta/611.
Full calving algorithm detailed [here](https://zenodo.org/records/10182710).