# Inefficient drainage system solver
## General Information
- **Solver Fortran File:** IDSSolver.f90
- **Solver Name:** IDSSolver
- **Required Output Variable(s):** IDSHead, IDSHead Residual, IDSHead Homologous and IDSHead Pressure
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** IDS Head Upper Limit (required if the upper limit is used)

## General Description
This solver treats the diffusion equation with a user-defined upper limit.

## SIF contents
The required keywords in the SIF file for this solver are given below. The IDSSolver can be used alone, coupling between the two layer is treated in the [EPLSolver](./HydrologyEPL.md) section

The hydrological system is only treated at the bed, it requires then a new body with a specific equation and intial condition, the Material and Body Force section are using the one from the ice.

```
Body 2
  Name = "hydrological system"
  Equation = 2
  Material = 1
  Body Force = 1
  Initial Condition = 2
End
! Initial condition for the hydrology
Initial Condition 2
  IDSHead = real 0.0
End
Constants
  Water Compressibility = Real 5.04e-4  !MPa-1
End
```
Only the parameters which are needed for the treatment of the hydrology are given here, you should add it to your existing Body Forces and Material.

```
Body Force 1
  Flow BodyForce 1 = 0.0
  Flow BodyForce 2 = 0.0
  Flow BodyForce 3 = -9.7696e15  ! or whichever value is used for gravity

  IDSHead Source Flux = Real 2.0  !water input into the sediment layer (distributed)
End
Material 1
! General Hydrology Parameters
  Water Density = Real MATC "1000.0*1.0E-06*(31556926.0)^(-2.0)"  !This is freshwater density

! IDS Solver
  IDS Transmitivity = Real 5.0e2
  IDS Porosity = Real 0.4
  IDS Thickness = Real 20.0
  IDS Compressibility = Real 1.0e-2

!Upper limit at the flotation limit
  IDSHead Upper Limit = Variable Depth, coordinate 3
      Real matc "tx(1)+tx(0)*0.91"
End
Solver 1
 Equation = "IDS Equation"

  Procedure = "ElmerIceSolvers" "IDSSolver"
  Variable = IDSHead
  Variable DOFs = 1

  Steady State Convergence Tolerance = Real 1.0E-5

  Linear System Solver = Direct
  Linear System Direct Method = umfpack
  Linear System Convergence Tolerance = Real 1.0E-7
  Linear System Residual Output = integer 1

  Nonlinear System Max Iterations = Integer 100
  Nonlinear System Convergence Tolerance = Real 1.0E-6
  Nonlinear System Relaxation Factor = Real 1.0

  !This deals with the upper limit (enabled if TRUE)
  Apply Dirichlet = Logical TRUE

  Exported Variable 1 = String "IDSHead Residual"
  Exported Variable 2 = String "IDSHead Homologous"
  Exported Variable 3 = String "IDSHead Pressure"

End
Equation 2
 Active Solvers (1) = 1
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

  IDSHead Flux BC = Logical True
  IDSHead Water Flux = Real 0.0
End

Boundary Condition 5
  Name = "Glacier snout"
  Target Boundaries = 5

! Take care to choose a value below or equals to the upper
1 limit of the water head

  IDSHead = variable coordinate 3, depth
    real matc "tx(0)+0.91*tx(1)"
End
```

## Example
Two basic tests can be found in [ELMER_TRUNK]/elmerice/Tests/Hydro_SedOnly and [ELMER_TRUNK]/elmerice/Tests/Hydro_Coupled.

## Reference
When used this solver can be cited using the following reference :
de Fleurian, B.; Gagliardini, O.; Zwinger, T.; Durand, G.; Le Meur, E.; Mair, D. & Råback, P. A double continuum hydrological model for glacier applications The Cryosphere, 2014, 8, 137-153
