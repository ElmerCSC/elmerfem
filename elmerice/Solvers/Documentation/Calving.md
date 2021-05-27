# Solver Calving2D

## General Information

- **Solver Fortran File:** Calving.F90
- **Solver Name:** Find_Calving
- **Required Output Variable(s):** Height
- **Required Input Variable(s):** List of Vars to interpolate (see example)
- **Optional Output Variable(s):** 
- **Optional Input Variable(s):** List of vars to nullify (set to zero)

## General Description
This is the solver for calving in 2D. It doesn't really solve any PDEs, but rather predicts if and where calving will occur, based on the 'crevasse depth calving law'. The crevasse depth calving law works by identifying connected groups of nodes where either surface or basal crevasses are predicted to be open based on the stress regime. If we are using the 'Basal Crevasse Model' (see SIF contents below), calving is predicted to occur where surface & basal crevasses meet. Otherwise, calving occurs where surface crevasses reach the waterline.

Calving location (x-coordinate) is computed, and then a vector defining the deformation of each node on the front is computed. In other words, displacing all the nodes on the front by the 'calving vector' produces the new calving front.

## SIF Contents
```
Solver 15
  Equation = Calving
  Procedure = "ElmerIceSolvers" "Find_Calving"
  Exec Solver = "After Timestep" !important

  Basal Crevasse Model = Logical True !Look for surface crevasses meeting waterline
  Old Way = Logical False

  Yield Stress = Real 0.05 !MPa (usually...)
  Basal FreeSurface = Logical True
  Basal FreeSurface Variable Name = String "Zs Bottom"

  ! the crevasse water depth parameters
  !Water Depth = Real $Dw
  ! Dw Start = Real $DwStart
  ! Dw Stop = Real $DwStop
  ! Dw Mode = String $DwMode

  Variable = Calving
  Variable DOFs = 2
  Exported Variable 1 = -dofs 1 "Calving Surface Index"
  Exported Variable 2 = -dofs 1 "Calving Basal Index"
  Exported Variable 3 = -dofs 1 "Crevasse Group ID"
End
```

# Solver Calving3D

## General Information
- **Solver Fortran File:** Calving3D.f90
- **Solver Name:** Find_Calving3D
- **Required Output Variable(s):** 
- **Required Input Variable(s):** 
- **Optional Output Variable(s):** 
- **Optional Input Variable(s):** 

## General Description

Information about the algorithm in this solver is provided [here](http://elmerfem.org/elmerice/wiki/doku.php?id=problems:calving).

In addition to the main solver, Calving3D makes use of two auxiliary solvers which must be present in the SIF: Calving Isosurface, ProjectCalving.

## SIF contents

```
Solver Options
Solver 19
  Equation = "3D Calving"
  Exec Solver = "After Timestep"
  Procedure = "ElmerIceSolvers" "Find_Calving3D"
  Solver Timing = Logical True

  Variable = String "Calving"
  Variable DOFs = 3

  Exported Variable 1 = -dofs 1 "CIndex"
  Exported Variable 1 DOFs = 1

  Calving Search Distance = Real 3000.0
  Calving Mesh Min LC = Real 30.0
  Calving Mesh Max LC = Real 100.0
  Calving Mesh LC Min Dist = Real 500.0
  Calving Mesh LC Max Dist = Real 1500.0

  Calving Append Name = String "$namerun"" !"
  Calving Move Mesh Dir = String "./results/"

  Project Calving Equation Name = String "CalvingProjection"
  Isosurface Equation Name = String "Calving Isosurface"
  Crevasse Penetration Threshold = Real 0.2 !this is the upper limit of average intact ice
  Minimum Calving Event Size = Real 1.0 !minimum front displacement length
  Pause Solvers Minimum Iceberg Volume = Real 1.0E6

  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations  = 2000
  Linear System Preconditioning = ILU1
  Linear System Convergence Tolerance = 1.0e-9
  Linear System Abort Not Converged = False
End

Solver 11
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

Solver 24
  Exec Solver = "Never" !auxiliary solver called by Calving3D
  Equation = "Calving Isosurface"
  Procedure = File "Isosurface" "IsosurfaceSolver"

  Isosurface Variable = String "ave_cindex"
  Isosurface Value = Real 0.0
End
```

## Examples
See [ELMER_TRUNK]/elmerice/Tests/Calving3D.
