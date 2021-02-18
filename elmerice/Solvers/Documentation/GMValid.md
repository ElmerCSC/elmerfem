# Solver GMValid
## General Information
- **Solver/User Function Fortran File:** GMValid.F90
- **Solver/User Function Name:** GMValid
- **Required Output Variable(s):** GMCheck
- **Required Input Variable(s):**
  - Grounded mask (GroundedMask Variable = String...) - this will default to 'GroundedMask' if not specified
- **Optional Output Variable(s):**
  - None
- **Optional Input Variable(s):**
  - None
- **Solver Keywords:**
  - None
  
## General Description
This is a very simple solver that’s more-or-less just a stripped-down copy of BasalMelt3D.F90 and exists to set up a mask variable for which ungrounded areas are connected to the fjord (i.e. the frontal boundary) and which aren’t.
This solver is part of the coupled calving-GlaDS-plumes suite described in [this document](./CoupledIceHydrologyCalvingPlumesDocumentation.md). It may require additional work to be used outside of this context.

## Known Bugs and Limitations
None

## SIF contents
{Provide examples of all the required SIF contents for a typical use of the solver/USF. Please ensure you include examples of lines that need to be added to other parts of the SIF (BCs, Material, etc.), not just the solver/USF block}
The required keywords in the SIF file for this solver are:

```
Solver 1
  Equation = "GMValid"
  Procedure = "ElmerIceSolvers" "GMValid"
  Variable = -dofs 1 "GMCheck"
  
  GroundedMask Variable = String "GroundedMask"
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
