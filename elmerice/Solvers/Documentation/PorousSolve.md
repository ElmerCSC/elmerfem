# Snow/firn rheology - Solver Porous Solver
## General Information
- **Solver Fortran File:** PorousSolve.f90
- **Solver Name:** PorousSolver
- **Required Output Variable(s):** Porous
- **Required Input Variable(s):** Relative Density
- **Optional Output Variable(s):** StrainRate, DeviatoricStress
- **Optional Input Variable(s):** None

## General Description
This solver computes the flow of snow/firn material (i.e., porous incompressible ice) using the snow/firn law proposed by Gagliardini and Meyssonnier (1997). The snow/ice rheological law is function of the relative density, which is a required input variable for this solver. The law depends on two function, a(D) and b(D), which are parametrized functions of the relative density D.

As optional output variables, strain-rate, deviatoric stress and spin can be calculated. In this solver, the nodal value of these quantities is computed as the average contribution from all the elements belonging at this given node. This can be slightly different than the evaluation obtained using the variational method (as in the [Strain-rate solver](./ComputeStrainRate.md) and the [ComputeDevStress](./ComputeDevStress.md) solver).

More details about the snow/firn law can be found here: [poroussolver.pdf](./poroussolver.pdf).

2024-01-23: The [ComputeDevStress](./ComputeDevStress.md) solver is now compatible with the Porous Solver. Using the latter is usually preferable than computing deviatoric stress directly in the Porous Solver. 

A New module PorousMaterialModels has been added in the Utils directory. It computes the effective viscosity and compressibility parameter to avoid code duplication (transparent for the user).

## SIF contents
The required keywords in the SIF file for the Porous Solver are:

```
! Define some useful parameters using MATC
$yearinsec = 365.25*24*60*60
$rhoi = 900.0/(1.0e6*yearinsec^2)
! B = 2 A, where A is the classical Glen's fluidity
$B = 20.0 ! MPa{^-3}a{^-1} T = -10°C
$n = 3.0
$gravity = -9.81*yearinsec^2

Constants
! give the name of the relative density variable 
  Density Name = String "Relative Density"
End

! this is the compressible Stokes solver
!----------------------------------------
Solver 1
  Equation = String "PorousFlow"
  Procedure = "ElmerIceSolvers" "PorousSolver"
  Variable = "Porous"
  Variable DOFs = 4 ! 4 in 3D (u,v,w,p) ; 3 in 2D (u,v,p)
  
  Optimize Bandwidth = False
! Use p elements
! Element = "p:1 b:4"
! Stablization Method = String pBubbles

  Exported Variable 1 = String "Relative Density"
  Exported variable 1 DOFs = Integer 1

! switch that in for post-processing issues only
   Exported Variable 2 = String "StrainRate"
   Exported variable 2 DOFs = Integer 6 ! 4 in 2D, 6 in 3D
   Exported Variable 3 = String "DeviatoricStress"
   Exported variable 3 DOFs = Integer 6 ! 4 in 2D, 6 in 3D
   Exported Variable 4 = String "Spin"
   Exported variable 4 DOFs = Integer 3 ! 1 in 2D, 3 in 3D

  Linear System Solver = 'Direct'
! Only Picard linearization available for this solver
  Nonlinear System Convergence Tolerance = 1.0E-05
  Nonlinear System Max Iterations = 50

  Steady State Convergence Tolerance = 1.0E-03
End

Equation 1
 Active Solvers(1) = 1
! Give the name of the porous solver variable
 Flow Solution Name = String "Porous"
End

! Gravity force is directly the ice density time the gravity
! It is further multiplied by the relative density in the Porous solver
Body Force 1
  Porous Force 1 = Real 0.0E00
  Porous Force 2 = Real 0.0E00
  Porous Force 3 = Real $gravity*rhoi 
End

Material 1
  Powerlaw Exponent = Real $n
  Min Second Invariant = Real 1.0E-10
  Fluidity Parameter = Real $B  ! MPa^{-3}a^{-1} 
 
! Just for output purpose, not needed by the Porous solver   
! Density as a function of relative density
  Density = Variable Relative Density
        Real MATC "tx*rhoi"
End

! Neumann type boundary condition
Boundary Condition 1
    Force 3 = Real -0.01
End

! or
Boundary Condition 1
    Normal Force = Real -0.01
End

! Dirichlet / Newton Boundary condition
! here: zero normal velocity and sliding
Boundary Condition 2
  Target Boundaries  = 2
  Normal-tangential Porous = True
  Porous 1 = Real 0.0
  Porous Slip Coeff 2 = Real 0.1
  Porous Slip Coeff 3 = Real 0.1
End
```

## Examples
An example using the Porous Solver can be found in [ELMER_TRUNK]/elmerice/examples/Test_Porous.

## Reference
The snow/firn rheological law is from:
Gagliardini O. and J. Meyssonnier, 1997. Flow simulation of a firn covered cold glacier. Annals of Glaciol., 24, p. 242-248.

Its implementation within Elmer/Ice and an application are presented in this reference:
Zwinger T. , R. Greve, O. Gagliardini , T. Shiraiwa and M. Lyly, 2007. A full Stokes-flow thermo-mechanical model for firn and ice applied to the Gorshkov crater glacier, Kamchatka. Annals of Glaciol., 45, p. 29-37.
