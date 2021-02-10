# User Function SourceCalcCalving
## General Information
- **User Function Fortran File:** USF_SourceCalcCalving.F90
- **User Function Name:** SourceCalc
- **Required Output Variable(s):** Source
- **Required Input Variable(s):**
  - SurfaceMelt (Surface Melt Variable Name = String...)
  - InternalMelt (Internal Melt Variable Name = String...)
- **Optional Output Variable(s):**
  - None
- **Optional Input Variable(s):**
  - None
- **Solver Keywords:**
  - Internal Melt = Logical... (Switch for whether you want to work out internal/basal melt)
  - Surface Melt = Logical... (Switch for whether you want to include surface melt)
  
## General Description
This is a USF that calculates the Hydraulic Potential Volume Source term required in the Body Force section of the SIF for GlaDS. If provided with a surface runoff variable (e.g. loaded in from a raster) and the temperature residual variable (for internal/basal melt), it will calculate the resulting internal melt and add on the surface melt for each node on the hydrology mesh (or, at the base of your 3D ice mesh, if you’re using GlaDS without all the other bells and whistles), so you can easily vary the source term spatially across your domain.
The USF also assumes the name of the hydraulic potential variable calculated by GlaDS is 'hydraulic potential', so don't change it.
This USF is part of the coupled calving-GlaDS-plumes suite described in [this document](../../Solvers/Documentation/CoupledIceHydrologyCalvingPlumesDocumentation.md). It may require additional work to be used outside of this context.

## Known Bugs and Limitations
None

## SIF contents
{Provide examples of all the required SIF contents for a typical use of the solver/USF. Please ensure you include examples of lines that need to be added to other parts of the SIF (BCs, Material, etc.), not just the solver/USF block}
The required keywords in the SIF file for this USF are:

```
Body Force 1
  Internal Melt = Logical True
  Surface Melt = Logical True

  !Set source to any old variable and then give actual variables below
  Hydraulic Potential Volume Source = Variable temp
    Real Procedure "ElmerIceUSFs" "SourceCalc"

  Internal Melt Variable Name = String "temp residual"
  Surface Melt Variable Name = String "runoff"
End
```
## Examples
TODO
An example in which the ... can be found here [ELMER_TRUNK]/elmerice/Tests/...

## References
None
