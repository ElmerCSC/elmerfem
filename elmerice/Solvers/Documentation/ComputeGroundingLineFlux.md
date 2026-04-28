# Solver ComputeGroundingLineFlux
## General Information
- **Solver  Fortran File:** ComputeGroundingLineFlux.F90
- **Solver/User Function Name:** ComputeGroundingLineFlux
- **Output Variable(s):**
    - tendligroundf: a global variable with the total flux accross the grounding line
    - ligroundf: a variable defined on the elements with the area mean local grounding flux 
- **Required Input Variable(s):**:
    - Flow Solution: the 3D velocity field used to compute the flux
    - GroundedMask: the grounded mask used to locate the gornding line
- **Optional Output Variable(s):** 
    - bc_area: the element area
- **Optional Input Variable(s):** None
- **Solver Keywords:** None

## General Description

Compute the ice flux accross the grounding line for a 3D simulation. The 3D velocity field will usually be solution of the Stokes equation. The groundind line is solution of the contact problem solved on the bottom boundary. The grounded mask defines the grounding line with groundedmask=0 the last grounded nodes.

This works only for vertically structured meshes. The grounding line flux is the flux through the vertical faces corresponding to edges where the two nodes have groundedmask=0. The variable *ligroundf* contains the area mean local grounding flux for each element that are partially grounded. The total flux is saved in the solver global variable *tendligroundf*.

The (projected) element area can be saved if the variable "bc_area" exists; tendligroundf=SUM_elements ligroundf*bc_area

The impelmentation is equivalent to the computation used in the SSA solver.

The solver must be executed on the bottom lowxer surface and the mesh must be vertically extruded.

As only be tested for wedge elements (extrusion of 2D triangles) but should work with bricks (verticall extrusion of 404 elements).

## SIF Contents

The required keywords in the SIF file for this solver are:

```
Solver ..
 ! Usually executed after time steps
  Exec Solver = After TimeStep
  Equation = "ComputeGlFlux"
  Procedure = "ElmerIceSolvers" "ComputeGroundingLineFlux"

  ! optionnal variable to save the (projected) element area
  Exported Variable 1 = -elem "bc_area"
End
```

## Examples
An example can be found in [ELMER_TRUNK]/elmerice/Tests/GroundingLineFlux

