# Hydrology - GlaDS model from Werder et al. (2013)
## General Information
- **Solver Fortran Files:** GlaDSCoupledSolver.F90 and GlaDSchannelSolver.F90
- **Solver Names:** GlaDSCoupledSolver, GlaDSsheetThickDummy and GlaDSchannelOut
- **Required Output Variable(s):** Hydraulic Potential, Sheet Thickness and Channel Area
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** Vclose, Wopen, Water Pressure, Effective Pressure, Sheet Discharge, Sheet Storage, Flux from Moulins (nodal variables) and Channel Flux (edge variable, to be exported in GlaDSchannelOut).
- **Optional Input Variable(s):** Zb

## Known bugs
If in parallel a moulin belongs to two partitions, the flux from the moulin is taken into account twice. This is may be not a bug, but the partitioning should avoid to have moulins at the boundary of partitions. The python tool makemoulin.py (in [ELMER_TRUNK]/elmerice/Meshers) take care to have no duplicated moulins on the boundaries of partitions.

If running the solver on a 3d internally extruded mesh, one should specify Preserve Edges = True and Preserve Baseline = True in the Simulation section.

If running on a “true” 3d mesh, the GlaDSchannelsSolver has to be executed on the 3d body (not the bottom surface body as for the two other GlaDS solvers).

## General Description
The complete description of the equations solved by the GlaDS solver can be found in Werder et al. (2013). The implementation follows exactly these equations, except that optionally the hydraulic potential can be computed at the top of the water sheet instead than at the bed (keyword: Neglect Sheet Thickness in Potential).

The GlaDS solver solves for the hydraulic potential, the water sheet thickness and the cross-sectional area of the channels. Whereas the two first variables are nodal variable and define continuous fields, the Channel area is a discrete field only defined on the edge of the elements.

The GlaDS model is composed of three solvers:
- GlaDSCoupledSolver is the main solver and couple the solve of the 3 main variables Hydraulic Potential, Sheet Thickness and Channel Area. Detail on the keywords for this solver are given below.
- GlaDSsheetThickDummy is just a solver to declare the Sheet Thickness variable as a primary variable.
- GlaDSchannelOut has two functions: declare that the Channel Area variable is an edge variable (Element = “n:0 e:1”) and create output vtu files for edge variables.
Since version Version 8.3 (Rev: b213b0c8), GlaDSchannelOut works for parallel simulation (no more vtk or acscii output, only vtu). These solvers only work in transient. They can be executed either on a 2d plane view mesh defining the bedrock or on the boundary of a 3d mesh. If using internal extrusion within Elmer see the [structured mesh](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:structuredmesh) page for essential keywords to preserve baseline and edges. More details about the specificity of the solvers are given below.

## SIF contents
The SIF examples given here are from the tests used in the [SHMIP](https://shmip.bitbucket.io/). The name of the variables as well as some constants have to be defined in the Constants section:

```
Constants
  Latent Heat = Real $Lw
  Gravity Norm = Real $gravity
  Fresh Water Density = Real $rhow
  Ice Density = Real $rhoi
  Sheet Thickness Variable Name = String "Sheet Thickness"
  Hydraulic Potential Variable Name = String "Hydraulic Potential"
  Channel Area Variable Name = String "Channel Area"
  Bedrock Variable Name = String "Zb"
End
```
The GlaDS solvers depend on a lot of physical parameters. The main parameters to be defined in the Material section are:

```
! For the sheet 
  Sheet Conductivity = Real $Ks
  Sheet flow exponent alpha = Real $alphas
  Sheet flow exponent beta = Real $betas
  Englacial Void Ratio = Real $ev
  Bedrock Bump Length = Real $lr
  Bedrock Bump High = Real $hr
  Sheet Closure Coefficient = Real $Ar
  
! For the Channels
  Channel Conductivity = Real $Kc
  Channel flow exponent alpha = Real $alphac
  Channel flow exponent beta = Real $betac
  Channel Closure Coefficient = Real $Ac
  Sheet Width Over Channel = Real $lc
  Pressure Melting Coefficient = Real $Ct
  Water Heat Capacity = Real $Cw

! Coupling with ice flow and glacier geometry
  Sliding Velocity = Real $ub
  Ice Normal Stress = Variable Coordinate 1
     Real MATC "rhoi*gravity*H(tx)"
```
In the Body Force section, one can set a water input source:

```
  Body Force 1
    Hydraulic Potential Volume Source = Real $Source
  End
```
GlaDSCoupledSolver solves for the three variables Hydraulic Potential, Sheet Thickness and Channel Area in a coupled way inside the solver itself (Coupled Max Iterations and Coupled Convergence Tolerance). Equations for the Hydraulic Potential and Channel Area are non linear. Only the equation for the Hydraulic Potential needs to solve a system. The two others are local and can be solved either explicitly, implcitely or using the Crank-Nicholson method.

```
Solver 1
  Equation = "GlaDS Coupled sheet"
  Procedure = "ElmerIceSolvers" "GlaDSCoupledSolver"
  Variable = -dofs 1 "Hydraulic Potential"

! activate or not the development of channels
  Activate Channels = Logical True           
  
! activate or not the growth of channels by melt            
  Activate Melt from Channels = Logical True 
  
! compute the hydraulic potential at the top of the water sheet (''False'') or at the bed (''True'')            
  Neglect sheet Thickness in Potential = Logical True  

! choices are EXPLICT, CRANK-NICOLSON, IMPLICIT
  Channels Integration method = String "Crank-Nicolson"
  Sheet Integration method = String "Implicit"

! define exported variables for visualization 
  Exported Variable 1 = -dofs 1 "Vclose"               ! closure velocity of the water sheet layer
  Exported Variable 2 = -dofs 1 "Wopen"                ! opening velocity of the water sheet layer     
  Exported Variable 3 = -dofs 1 "Water Pressure"       ! water pressure at the base 
  Exported Variable 4 = -dofs 1 "Effective Pressure"   ! effective pressure at the base
  Exported Variable 5 = -dofs 2 "Sheet Discharge"      ! water discharge (vector) in the water sheet layer 
  Exported Variable 6 = -dofs 1 "Sheet Storage"        ! storage in the water sheet layer
  Exported Variable 8 = -dofs 1 "Zs"                   
  Exported Variable 9 = -dofs 1 "Zb"
  Exported Variable 10 = -dofs 1 "Flux from Moulins"   ! flux of water from the moulins (Qm)  

  Linear System Solver = Direct
  Linear System Direct Method = umfpack

  Nonlinear System Max Iterations = 30
  Nonlinear System Convergence Tolerance  = 1.0e-6
  Nonlinear System Relaxation Factor = 1.00

  Coupled Max Iterations = Integer 10
  Coupled Convergence Tolerance = Real 1.0e-3

  Steady State Convergence Tolerance = 1.0e-03
End
```
GlaDSsheetThickDummy is just here to declare the variable Sheet Thickness

```
Solver 2
  Equation = "GlaDS Thickness sheet"
  Procedure = "ElmerIceSolvers" "GlaDSsheetThickDummy"
  Variable = -dofs 1 "Sheet Thickness"
End
```
GlaDSchannelOut allows to declare the Channel Area variable as an edge variable and to save edge variables. Output file formats is VTU.

```
Solver 3
  Exec Solver = After Saving
  Equation = "GlaDS Channel OutPut"
  Procedure = "ElmerIceSolvers" "GlaDSchannelOut"
  Variable = -dofs 1 "Channel Area"
! the variable is define on the edges only
  Element = "n:0 e:1"

  Exported Variable 1 = -dofs 1 "Channel Flux"

  VTU OutPutFile = Logical False   ! set to TRUE to have VTU output
  VTU BinaryFile = Logical False   

  Channels OutPut Directory Name = String "results"
  Channels OutPut File Name = String "$namerun"_channels"
End
```
The possible boundary conditions are:
- channels not allowed to grow (recommended on all the domain boundary)
```
Boundary Condition 1
  Target Boundaries(2) = 1 3
  No Channel BC = Logical True
End
```
- fixed value of the Hydraulic Potential at some specific outlet nodes:
```
Boundary Condition 2
  Name = "point front"
  Target Coordinates(1,2) = Real 0.0 0.0
  Hydraulic Potential = Variable Zb
     Real MATC "rhow*gravity*tx"
End
```
- the possibility to impose water flux at some nodes in the domain (moulins type inflow). The nodes have to be declared as node element (101) in the mesh (mesh.header and mesh.boundary files have to be modified by hand for that, or using the python tool makemoulin.py - see below).
```
Boundary Condition 4
  Name = "moulins"
  Target Boundaries(20) = 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
  Moulin Storage = Logical True
  Moulin Area = Real $Am
  Moulin Flux = Real $4.5*yearinsec
End
```
## Making a mesh with Moulins
Moulins are 101 boundary elements. ElmerGrid does not export correctly 101 boundary elements from gmsh or when partitioning a mesh. To add 101 boundary elements to an existing mesh, thanks to Mondher Chekki (IGE), one can use the python tool makemoulin.py (in [ELMER_TRUNK]/elmerice/Meshers).

Usage:

`python makemoulin.py --meshdir mesh_dir --moulin moulin_file --partition number_of_partition`
where moulin_file is an ascii file which contains the (x,y) coordinates of the moulins. The same file has to be used in gmsh so that nodes with the moulin coordinates already exist.

## Examples
Examples using the GlaDS Solver can be found in [ELMER_TRUNK]/elmerice/Tests/.

## References
The description of the GlaDS model is in:
- Werder M.A., I.J. Hewitt, C.G. Schoof and G.E. Flowers, 2013. Modeling channelized and distributed subglacial drainage in two dimensions. Journal of Geophysical Research: Earth Surface, 118(4), 2140-2158.
The implementation of the GlaDS model in Elmer/Ice is described here:
- Gagliardini O. and M. Werder, 2018. Influence of increasing surface melt over decadal timescales on land-terminating Greenland-type outlet glaciers, Journal of Glaciology, 64(247), 700-710, doi:10.1017/jog.2018.59
Results of the SHMIP experiments using Elmer/Ice are discussed in the SHMIP paper:
- De Fleurian, B., M. Werder, S. Beyer, D. Brinkerhoff, I. Delaney, C. Dow, C., J. Dows, O. Gagliardini, M.J. Hoffman, R. LeB Hooke, J. Seguinot, A.N. Sommers, 2018. SHMIP The subglacial hydrology model intercomparison Project. Journal of Glaciology, 1-20. doi:10.1017/jog.2018.78
