# Solver BasalMelt3D
## General Information
- **Solver Fortran File:** BasalMelt3D.F90
- **Solver/User Function Name:** BasalMelt3D
- **Required Output Variable(s):** BasalMeltRate
- **Required Input Variable(s):** 
  - Grounded mask (GroundedMask Variable = String...) - this will default to 'GroundedMask' if not specified
- **Optional Output Variable(s):** 
  - None
- **Optional Input Variable(s):**
  - None
- **Solver Keywords:**
  - Basal Melt Stats File = String... (The path to write a file containing some stats to)
  - Basal Melt Mode = String... (can be 'seasonal' or 'off'. The latter is fairly self-explanatory; the former requires the following additional options:)
  - Basal Melt Summer Rate = Real... (The rate to apply basal melt at in summer)
  - Basal Melt Winter Rate = Real... (The rate to apply basal melt at in winter)
  - Basal Melt Summer Start = Real... (The time in the simulation to begin using summer melt rates (expressed as a number between 0 and 1))
  - Basal Melt Summer Stop = Real... (The time in the simulation to stop using summer melt rates (expressed as a number between 0 and 1))
  
  
## General Description
This is a very simple solver written by Joe Todd that applies a specified basal melt rate to any ungrounded parts of the glacier base.
This solver is part of the coupled calving-GlaDS-plumes suite described in [this document](./CoupledIceHydrologyCalvingPlumesDocumentation.md). It may require additional work to be used outside of this context.

## Known Bugs and Limitations
None

## SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 1
  Equation = "Basal Melt"
  Procedure = File "ElmerIceSolvers" "BasalMelt3D"

  Variable = "BasalMeltRate"
  Basal Melt Stats File = String "basalmeltstats.txt"
  GroundedMask Variable = String "GroundedMask"

  Basal Melt Mode = String $basalmeltmode
  Basal Melt Summer Rate = Real 100.0
  Basal Melt Winter Rate = Real 10.0
  Basal Melt Summer Start = Real 0.4
  Basal Melt Summer Stop = Real 0.6
End

Boundary Condition 1
  Calving Front Mask = Logical True
End
```
## Examples
TODO
An example in which the ... can be found here [ELMER_TRUNK]/elmerice/Tests/...

## References
None
