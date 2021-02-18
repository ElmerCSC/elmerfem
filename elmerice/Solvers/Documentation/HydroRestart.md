# Solver HydroRestart
## General Information
- **Solver Fortran File:** HydroRestart.F90
- **Solver/User Function Name:** HydroRestart
- **Required Output Variable(s):** None
- **Required Input Variable(s):**
  - None
- **Optional Output Variable(s):** 
  - None
- **Optional Input Variable(s):**
  - {e.g. SomeInputVar (Extra Input Name = String...)}
  - hp: Restart Variable 1 = String... (any GlaDS variables that belong with the main GlaDS solver subroutine (GlaDSCoupledsolver in GlaDSCoupledSolver.F90))
  - channel: Restart Variable 1 = String... (any GlaDS channel variables)
  - sheet: Restart Variable 1 = String... (any GlaDS variables that belong with the dummy sheet thickness solver (GlaDSsheetThickDummy in GlaDSCoupledSolver.F90); in practice, this is just sheet thickness)
- **Solver Keywords:**
  - None
  
## General Description
Largely a direct copy of the Restart() subroutine within the Elmer source code, with a few modifications to disentangle it from all the other restart machinery and to make it pick up the right variables from the right place and send them to the right place on the secondary hydrology mesh. Assumes the hydraulic potential, channel area and sheet thickness variables are all called by those default names, so don't change them. Uses the namespace feature in Elmer.
This solver is part of the coupled calving-GlaDS-plumes suite described in [this document](./CoupledIceHydrologyCalvingPlumesDocumentation.md). It may require additional work to be used outside of this context.

## Known Bugs and Limitations
None

## SIF contents
The required keywords in the SIF file for this solver/USF are:

```
Solver 1
  Equation = "HydroRestart"
  Procedure = File "ElmerIceSolvers" "HydroRestart"
  Exec Solver = "Before Simulation"

  hp: Restart Variable 1 = String "Hydraulic Potential"
  hp: Restart Variable 2 = String "Effective Pressure"
  hp: Restart Variable 3 = String "Water Pressure"
  hp: Restart Variable 4 = String "Sheet Discharge"
  hp: Restart Variable 5 = String "Sheet Storage"
  hp: Restart Variable 6 = String "wopen"
  hp: Restart Variable 7 = String "vclose"

  channel: Restart Variable 1 = String "channel flux"
  channel: Restart Variable 2 = String "channel area"

  sheet: Restart Variable 1 = String "sheet thickness"
End

```
## Examples
TODO
An example in which the ... can be found here [ELMER_TRUNK]/elmerice/Tests/...

## References
None
