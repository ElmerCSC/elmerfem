# Solver Flotation
## General Information
- **Solver Fortran File:** Flotation.F90
- **Solver Name:** Flotation
- **Required Output Variable(s):** 
	- Zb (or Variable name prescribed by *Bottom Surface Name*)
	- Zs (or Variable name prescribed by *Top Surface Name*)
- **Required Input Variable(s):** H (or Variable name prescribed by *Thickness Variable Name*)
- **Optional Output Variable(s):** 
	- GroundedMask
	- sftgif,sftgrf,sftflf
- **Optional Input Variable(s):** bedrock

## History
- Rev. 698263163: add computation of area fractions
	- New Solver Keyword " compute ice area fractions = Logical"

## General Description

### Description 

The aim of this solver is to apply the flotation criterion to compute the top and bottom surface elevation, knowing the ice thickness. In general it will be used with the [SSA Solver](./SSA.md) and [Thickness Solver](./ThicknessSolver.md)

- The bottom surface elevation $z_b$ is computed as:

  $z_b=z_{sea} - H \rho_i / \rho_w$  
  where  
  - $H$ is the ice thickness,   
  - $z_{sea}$ is the (sea) water level elevation,   
  - $\rho_i$ is the mean ice density   
  - $\rho_w$ is the (sea) water density.

- If the bedrock variable is present, $z_b=max(z_b,bedrock)$  

- The top surface elevation $z_s$ is then simply given as $z_s = z_b + H$

- If the GroundedMask variable is present:  
  - GroundedMask=1 where $z_b=bedrock$ (grounded ice)
  - GroundedMask=-1 where $z_b>bedrock$ (floating ice)
  - GroundedMask=0 at the grounding line (list of nodes where $z_b=bedrock$ but the nodes belong to at least one grounded (all nodes grounded) and one floating (at leat one node floating) element  

- IF compute ice area fractions = TRUE, compute the element area fractions sftgif (land_ice_area_fraction), sftgrf (grounded_ice_sheet_area_fraction), sftflf (floating_ice_shelf_area_fraction). No sub-scheme is used so values are 0._dp or 1._dp. An element is considered floating if at least one node is floating, otherwise it is grounded. If the solution for the Thickness is limited, then if all nodal H == Lower Limit, all teh area fractions are set to 0.

### Remarks  

- It might be interesting to compute the top and bottom surface elevation rate of change; this can be done using internal Elmer functionality with the keyword,e.g. *Zs Calculate Velocity = Logical True* in the solver where *Zs* is created as an exported variable. See [ElmerSolver Manual](http://www.nic.funet.fi/pub/sci/physics/elmer/doc/ElmerSolverManual.pdf), section **13.4 Exported and derived variables**

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

  ![OPTIONAL :] 
  Bottom Surface Name = String "zb" ![Default: Zb]  
  Top Surface Name = String "zs" ![Default: Zs]  
  Thickness Variable Name = String "H" ![Default: H] 

   Exported Variable 1 = -dofs 1 "Zs"
   Exported Variable 2 = -dofs 1 "Zb"
   Exported Variable 3 = -dofs 1 "bedrock"

  ![OPTIONAL :] rates of changes of Zs and zb can be computed with 
  ! (in the solver where they are created as Exported Variables):
  ! Zs Calculate Velocity = Logical True
  ! Zb Calculate Velocity = Logical True

  ![OPTIONAL :]
  compute ice area fractions = Logical TRUE
  
  Exported Variable 1 = -elem "sftgif"
  Exported Variable 2 = -elem "sftgrf"
  Exported Variable 3 = -elem "sftflf"
End
```

## Examples
For examples look in your elmer source distribution under
[ELMER_TRUNK]/elmerice/Tests/SSA_IceSheet
