# Solver CalvingHydroInterp
## General Information
- **Solver Fortran File:** CalvingHydroInterp.F90
- **Solver Name:** IceToHydroInterp and HydroToIceInterp
- **Required Output Variable(s):** None
- **Required Input Variable(s):**
  - Normal stress (assumed to be called 'normalstress')
  - Velocity (assumed to be called 'velocity')
  - Temperature residual (assumed to be called 'temp residual')
  - Water pressure (assumed to be called 'water pressure')
  - Effective pressure (assumed to be called 'effective pressure')
- **Optional Output Variable(s):**
  - None
- **Optional Input Variable(s):**
  - Grounded mask (assumed to be called 'GroundedMask')
  - Grounded mask validity mask (assumed to be called 'GMCheck')
- **Solver Keywords:**
  - (All of these are for IceToHydroInterp; there are no options for the other subroutines)
  - Reference Node (3) = Integer x y z (The reference node, usually somewhere on the front, that you want to use as a distance marker in artefact correction)
  - Threshold Distance = Real... (The distance from the reference node beyond which you want to force groundedness (removes GroundedMask artefacts inland))
  - Side = Logical... (This isn’t a solver option, but something that should be set in the boundary condition section of the SIF. The interpolation routine tends to create a lot of artefacts on the lateral boundaries of the domain – if you set `Side = Logical True` in the boundary condition sections for the sidewalls of the hydrology mesh, this will force them to be grounded and remove the frequent ungrounded artefacts.)
  
## General Description
This solver handles the interpolation of necessary variables between the ice and hydrology meshes. It also corrects interpolation artefacts that will otherwise nix your simulation sooner or later, and ensures conservation of the temperature residual (one of the interpolated variables) to stop the glacier accidentally destroying or creating some energy….
The file contains two main subroutines, imaginatively titled “IceToHydroInterp” and “HydroToIceInterp”. Make sure you get them the right way round. IceToHydroInterp interpolates the ice normal stress, velocity, grounded mask and temperature residual over to the hydrology mesh and then spends a lot of time clearing up artefacts and conserving the temperature residual. HydroToIceInterp is much simpler, as the hydrology mesh is usually finer than the ice mesh, so the interpolation routine doesn’t create anywhere near as many artefacts in problematic locations. Therefore, it pretty much just interpolates the water pressure and effective pressure onto the ice mesh.
There are also two small subroutines: “HydroWeightsSolver” and “IceWeightsSolver”. These calculate the boundary weights used in the main routines (the reason this happens in a separate solver is complicated – suffice to say it does exist). These need to be called as solvers before the relevant interpolation routines (IceToHydro or HydroToIce) are called, otherwise they’ll crash. IceWeightsSolver needs to run every time the ice mesh is updated (probably every n timesteps); HydroWeightsSolver every time the hydrology mesh is updated (probably never, so it can just run once at the start of the simulation).
This solver is part of the coupled calving-GlaDS-plumes suite described in [this document](./CoupledIceHydrologyCalvingPlumesDocumentation.md). It may require additional work to be used outside of this context.

## Known Bugs and Limitations
None

## SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 1
 Equation = "IceToHydroInterp"
  Procedure = "ElmerIceSolvers" "IceToHydroInterp"

  Load Reader Variables = Logical True
  Number Of Variables To Read = Integer 3
  Reader Solver 1 = Integer 4
  Reader V1 = String "zb"
  Reader Solver 2 = Integer 5
  Reader V2 = String "runoff"
  Reader Solver 3 = Integer 53
  Reader V3 = String "hydroweights"

  Reference Node(3) = Integer -209278 -2135324 0
  Threshold Distance = Real 10000.0
End

Solver 2
  Equation = "HydroToIceInterp"
  Procedure = "ElmerIceSolvers" "HydroToIceInterp"
End

Solver 3
  Equation = IceWeightsCalculator
  Procedure = "ElmerIceSolvers" "IceWeightsSolver"
  Variable = -dofs 1 "IceWeights"
End

Solver 4
  Equation = HydroWeightsCalculator
  Exec Solver = "Before Simulation"
  Procedure = "ElmerIceSolvers" "HydroWeightsSolver"
  Mesh = "." "mesh/HydroMeshFinalNO"
  Variable = -dofs 1 "HydroWeights"
End

Boundary Condition 1
  Side = Logical True !if a lateral boundary; set to False if front or inflow
End
```
## Examples
TO DO
An example in which the ... can be found here [ELMER_TRUNK]/elmerice/Tests/...

## References
None
