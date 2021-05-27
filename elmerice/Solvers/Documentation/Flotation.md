# Solver Flotation
## General Information
- **Solver Fortran File:** Flotation.F90
- **Solver Name:** Flotation
- **Required Output Variable(s):** Zb (or Variable name prescribed by Bottom Surface Name), Zs (or Variable name prescribed by Top Surface Name)
- **Required Input Variable(s):** H (or Variable name prescribed by Thickness Variable Name)
- **Optional Output Variable(s):** GroundedMask, bedrock, DZbDt, DZsDt
- **Optional Input Variable(s):** None

## General Description
The aim of this solver is to apply the flotation criterion to compute the top and bottom surface elevation, knowing the ice thickness. In general it will be used with the [SSA Solver](./SSA.md) and [Thickness Solver](./ThicknessSolver.md)

The bottom surface elevation z_b is computed as:

*z_b=z_sea - H rho_i / rho_w* where *H* is the ice thickness, *z_sea* is the (sea) water level elevation, *rho_i* is the mean ice density and *rho_w* is the (sea) water density.

- If the bedrock variable is present, *z_b=max(z_b,bedrock)*
- If the GroundedMask variable is present:
  - GroundedMask=1 where *z_b=bedrock* (grounded ice)
  - GroundedMask=-1 where *z_b>bedrock* (floating ice)
  - GroundedMask=0 at the grounding line (list of nodes where *z_b=bedrock* but the nodes belong to at least one grounded (all nodes grounded) and one floating (at leat one node floating) element
The top surface elevation *z_s* is then simply given as *z_s = z_b + H*

If the variable DZbDt and/or DZsDt are present then the elevation change is given as *dz/dt={Delta z} / {Delta t}*

## SIF contents
```
Constants
 Sea level = Real .... !z_sea
 Water Density = Real ... !rho_w 
End

Material 1
  SSA Mean Density = Real  ... !rho_i 
End

Solver 3
   Equation = "Flotation"
   Procedure = "ElmerIceSolvers" "Flotation"
   Variable =  "GroundedMask"

   Exported Variable 1 = -dofs 1 "Zs"
   Exported Variable 2 = -dofs 1 "Zb"
   Exported Variable 3 = -dofs 1 "bedrock"
End
```

## Examples
For examples look in your elmer source distribution under
[ELMER_TRUNK]/elmerice/Tests/SSA_IceSheet
