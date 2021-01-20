# Efficient drainage system solver
## General Information
- **Solver Fortran File:** EPLSolver.f90
- **Solver Name:** EPLSolver
- **Required Output Variable(s):** EPLHead, Open EPL and EPLHead Homologous
- **Required Input Variable(s):** IDSHead Residual and IDSHead
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** EPL Head Upper Limit (required if the upper limit is used)

## General Description
This solver treat the diffusion equation for the efficient drainage system. The activation of the layer is dealt by the Open EPL mask. The upper limit is needed to increase the size of the drainage system once the water head in the efficient layer reaches a given value.

## SIF contents
The required keywords in the SIF file for this solver are given below. The [IDSSolver](./HydrologyIDS.md) is needed to use the efficient layer solver. In the following, we assume that the SIF file is written as shown in the [IDSSolver](./HydrologyIDS.md) section; the SIF should then be written as follows.

```
! Initial condition for the hydrology
Initial Condition 2
  EPLHead = real 0.0
End
```
The transfer between the two layers is dealt with here. An example is given here with the [WaterTransfer](../../UserFunctions/Documentation/WaterTransfer.md) user function, but one could choose to define an other type of transfer (just keep in mind the sign of the transfer, positive if the flux is from the efficient to the inefficient system).

```
Body Force 1
EPLToIDS Transfer = Variable Coordinate 1
  Real Procedure "ElmerIceUSF" "EPLToIDS"
    
EPLHead Passive = Equals Open EPL
Active Element Min Nodes = Integer 1
End
Material 1

!this is a variable of the transfer USF
  Leakage Factor = Real 20.0

! EPL Solver
  EPL Transmitivity = Real 2.5e5
  EPL Porosity = Real 0.4
  EPL Thickness = Real 1.0
  EPL Compressibility = Real 1.0e-2

!Upper limit at the flotation limit
  EPLHead Upper Limit = Variable Depth, coordinate 3
      Real matc "tx(1)+tx(0)*0.91"
End
Solver 2
 Equation = "EPL Equation"

  Procedure = "ElmerIceSolvers" "EPLSolver"
  Variable = EPLHead
  Variable DOFs = 1

  Steady State Convergence Tolerance = Real 1.0E-5

  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Linear System Convergence Tolerance = Real 1.0E-7
  Linear System Residual Output = integer 1

  Nonlinear System Max Iterations = Integer 100
  Nonlinear System Convergence Tolerance = Real 1.0E-6
  Nonlinear System Relaxation Factor = Real 1.0

  IDS Residual Name = String "IDSHead Residual"
  IDS Load Name = String "IDSHead"

  !This deals with the upper limit (enabled if TRUE)
  Apply Dirichlet = Logical TRUE

  Exported Variable 2 = String "EPLHead Homologous"
  Exported Variable 3 = String "Open EPL"

End
Equation 2
 Active Solvers (2) = 1 2
End
```
The boundary condition of the hydrological model should be applied on a 1D boundary located at the corner between the side and bed of the mesh.

```
Boundary Condition 4
  Name = "Lower frame"
  Target Boundaries = 4
  
! Flux condition on the borders of the hydrological domain
! Zero flux is not a necessary input as it is the natural 
! boundary condition of the system

  EPLHead Flux BC = Logical True 
  EPLHead Water Flux = Real 0.0
End

Boundary Condition 5
  Name = "Glacier snout"
  Target Boundaries = 5
  
! Take care to choose a value bellow or equals to the upper 
1 limit of the water head
  
  EPLHead = variable coordinate 3, depth
    real matc "tx(0)+0.91*tx(1)"
End
```

## Example
Two basic tests can be found in [ELMER_TRUNK]/elmerice/Tests/Hydro_SedOnly and [ELMER_TRUNK]/elmerice/Tests/Hydro_Coupled.

## Reference
When used this solver can be cited using the following reference :
de Fleurian, B.; Gagliardini, O.; Zwinger, T.; Durand, G.; Le Meur, E.; Mair, D. & Råback, P. A double continuum hydrological model for glacier applications The Cryosphere, 2014, 8, 137-153
