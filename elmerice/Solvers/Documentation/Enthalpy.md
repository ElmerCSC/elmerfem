# Enthalpy Solver
## General Information
- **Solver Fortran File:** EnthalpySolver.f90
- **Solver Name:** EnthalpySolver
- **Required Output Variable(s):** Enthalpy_h, Phase Change Enthalpy, Temperature and Water Content
- **Required Input Variable(s):** a velocity field
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
Solves the enthalpy equation:

*rho {{\partial H}/{\partial t}} + rho u .grad H = div(k grad H) + tr (sigma epsilon) + Q_lat*

where

- *H* is the enthalpy variable
- *rho* is the ice density
- *u* is the ice velocity vector
- *k* is the enthalpy diffusivity
- *tr (sigma epsilon)* the strain heating
- *Q_lat* a complementary source term accounting for melt water refreezing
Enthalpy is defined as a function of the water content omega and the temperature T, such that:

*If H < H_f, then H(T, omega) = int_T0^T C_p (T) dT*

*If H > H_f, then H(T, omega) = int_T0^Tm C_p (T) dT + omega L*

where

- *H_f* is the enthalpy of fusion, defined from the fusion temperature according to the pressure dependent Clausius-Clapeyron relationship.
- *C_p* is the temperature dependant heat capacity, defined as C_p = AT+B
- *L* is the latent heat of fusion
For the boundary conditions, a flux (Enthalpy Heat Flux) has the same meaning than for the temperature solver (W/m2). For a Dirichlet boundary condition on the enthalpy variable, the same definition as in the solver has to be used, i.e. H(T, omega) = int_T0^T C_p (T) dT. See example below.

## SIF contents
In this example, ice velocity is in m/s and pressure in MPa.

```
Solver 2
  Equation = String "Enthalpy Equation"
  Procedure = File "ElmerIceSolvers" "EnthalpySolver"
  Variable = String "Enthalpy_h"
  Linear System Solver = "Iterative"
  Linear System Iterative Method = "BiCGStab"
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0E-07
  Linear System Abort Not Converged = True
  Linear System Preconditioning = "ILU0"
  Linear System Residual Output = 1
  Steady State Convergence Tolerance = 1.0E-04
  Nonlinear System Convergence Tolerance = 1.0E-03
  Nonlinear System Max Iterations = 10
  Nonlinear System Relaxation Factor = Real 1.0

  Apply Dirichlet = Logical True
  Stabilize = True

  Exported Variable 1 = String "Phase Change Enthalpy" ! (J kg-1)
  Exported Variable 1 DOFs = 1

  Exported Variable 2 = String "Water Content" ! (%)
  Exported Variable 2 DOFs = 1

  Exported Variable 3 = String "temperature" ! (°C)
  Exported Variable 3 DOFs = 1
End

Constants
 T_ref_enthalpy = real 200.0 !(J kg-1)
 L_heat = real 334000.0 !(J kg-1)
 ! Cp(T) = A*T + B
 Enthalpy Heat Capacity A = real 7.253 !(J kg-1 K-2)
 Enthalpy Heat Capacity B = real 146.3 !(J kg-1 K-1)
 P_triple = real 0.061173 !Triple point pressure for water (MPa)
 P_surf = real 0.1013 ! Surface atmospheric pressure(MPa)
 beta_clapeyron = real 0.0974 ! clausus clapeyron relationship (K MPa-1)
End

Body Force 1
  Heat Source = real 0.0
End

Material 1
  Enthalpy Density = real 917.0 !(kg m-3)
  Enthalpy Heat Diffusivity = Real $2.1/2050.0 ! = k / Cp (kg m-1 s-1)
  Enthalpy Water Diffusivity = real 1.045e-4 ! (kg m-1 s-1)
End

! bed rock interface
Boundary Condition 1
  Target Boundaries = 1
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
  Velocity 3 = Real 0.0

  Enthalpy Heat Flux BC = logical True
  Enthalpy Heat Flux = real 0.02 !(W m-2)
End

! Upper Surface
Boundary Condition 2
  Target Boundaries = 2
  Enthalpy_h = variable coordinate 3
    real MATC "25000.0/150.0*(tx-3250)+140000.0" ! (J kg-1)
End
End
```

## Examples
An example solving for the enthalpy within the Tete Rousse glacier assuming an elevation dependent enthalpy at the upper surface can be found in [ELMER_TRUNK]/elmerice/Tests/Enthalpy.

## References
Gilbert, A., O. Gagliardini, C. Vincent, and P. Wagnon, 2014. A 3-D thermal regime model suitable for cold accumulation zones of polythermal mountain glaciers, J. Geophys. Res. Earth Surf., 119, doi:10.1002/2014JF003199.
