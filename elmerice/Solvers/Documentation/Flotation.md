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
    - Haf,Haf0
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
  - GroundedMask=0 at the grounding line (list of nodes where $z_b=bedrock$ but the nodes belong to at least one grounded (all nodes grounded) and one floating (at least one node floating) element  

- If the **Haf** variable is present, compute the nodal heigh above flotation as:  
    - $Haf=H$ if $bedrock > zsea$
    - $Haf=H - (zsea-bedrock)*rhow/rhoi$ if grounded (GroundedMask=[0,1])
    - $Haf=0$ if floating (GroundedMask=-1)

- If the **Haf0** variable is present :
    - $Haf0=H$ if $bedrock > zsea$
    - $Haf0=H - (zsea-bedrock)*rhow/rhoi$
    - => the 0-isocontour of Haf0 exactly gives the sub-element grounding line location, and is used by the SEP2 sub-element integration for the ice fraction areas (see below)


- IF compute ice area fractions = TRUE, compute the element area fractions :
    - sftgif (land_ice_area_fraction)
    - sftgrf (grounded_ice_sheet_area_fraction)
    - sftflf (floating_ice_shelf_area_fraction)
    - sftgif=sftgrf+sftflf
    - If the solution for the Thickness is limited, then if all nodal H == Lower Limit, all the area fractions are set to 0.
    - To compute sftgrf and sftflf the flotation criterion is directly evaluated at the IPs:
        - the number of IPs at partly grounded elements can be changed with the solver kw 'GL integration points number = Integer N'. This is SEP3 in Seroussi et al. (2014). See [ELMER_TRUNK]/fem/src/Integration.F90 for the possible integration rules.
        - Adaptive integration where the elements are splitted following the sub-element GL location is activated using the solver kw 'Sub-Element GL parameterization = logical True'. This is SEP2 in Seroussi et al. (2014) and will give the exact grounded (resp. floating) area. The variable **Haf0** is required to locate the GL within the elements. 

Seroussi, H., Morlighem, M., Larour, E., Rignot, E., and Khazendar, A.: Hydrostatic grounding line parameterization in ice sheet models, The Cryosphere, 8, 2075–2087, https://doi.org/10.5194/tc-8-2075-2014, 2014.

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

  ![OPTIONAL :] Height above flotation
   Exported Variable 4 = Haf
   Exported Variable 5 = Haf0
  

  ![OPTIONAL :] rates of changes of Zs and zb can be computed with 
  ! (in the solver where they are created as Exported Variables):
  ! Zs Calculate Velocity = Logical True
  ! Zb Calculate Velocity = Logical True

  ![OPTIONAL :]
  compute ice area fractions = Logical TRUE

  ! activate Sub element integration:
  ! SEP2
  Sub-Element GL parameterization = Logical True 
  !SEP3
  ! 20 is the maximum coded inetgration rule for triangles
  ! GL integration points number = Integer 20
  
  Exported Variable 6 = -elem "sftgif"
  Exported Variable 7 = -elem "sftgrf"
  Exported Variable 8 = -elem "sftflf"
End
```

## Examples
For examples look in your elmer source distribution under :
- [ELMER_TRUNK]/elmerice/Tests/SSA_IceSheet  
- [ELMER_TRUNK]/elmerice/Tests/Flotation 
