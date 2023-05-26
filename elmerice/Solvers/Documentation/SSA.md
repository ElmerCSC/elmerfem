# Solver Shallow Stream/Shelf Approximation (SSA)
## General Information
**:!: Important changes have been made in SSABasalSolver. This doc applies from Rev. 6480.**

- **Solver Fortran File:** SSASolver.f90
- **Solver Name:** (1) SSABasalSolver, (2) GetMeanValueSolver and (3) SSASolver
- **Required Output Variable(s):**
  - (1) SSAVelocity
  - (2) Mean Viscosity and Mean Density
  - (3) SSAFlow
- **Required Input Variable(s):**
  - (1) Zb, Zs and Effective Pressure when using the Coulomb type friction law
  - (2) Depth
  - (3) Depth, FreeSurfGrad1, FreeSurfGrad2 and SSABasalFlow
- **Optional Output Variable(s):** 
  - (1) strbasemag : element-average basal friction
  - (2) Ceff: nodal effective friction coefficient
- **Optional Input Variable(s):** None

## History
- Rev. 26a079785: introduce "regularised Coulomb" friction law
	- Update Keyword: MATERIAL : SSA Friction Law = String "regularised Coulomb"
	- New Keyword: MATERIAL : SSA Friction Threshold Velocity = Real

- Rev. afcafe865: move friction law in a separate module
	- New module SSAMaterialModels in the Utils dir. with the friction laws

- Rev. afcafe865:  
	- Add post-processing options:  
		- compute element-average basal friction if variable "strbasemag" is found and is by element  
		- compute nodal effective friction coefficient if variable "Ceff" is found

## General Description
### Ice flow
The SSABasalSolver solves the classical SSA equation, it has been modified in Rev. 6440 to be executed either on a grid of dimension lower than the problem dimension itself (i.e. the top or bottom grid of a 2D or 3D mesh for a SSA 1D or 2D problem), or on a grid of the same dimension as the problem (i.e. 2D mesh for a 2D plane view SSA solution).
It will work on a 3D mesh only if the mesh has been extruded along the vertical direction and if the baseline boundary conditions have been preserved (to impose neumann conditions).
The mandatory input variables are the bottom surface elevation and top surface elevation variables called Zb and Zs, respectively.
For the Flow law the SSA solver uses a “power-law” formulation and uses the keywords Viscosity Exponent, Critical Shear Rate, and Mean Viscosity. It doesn't work with the built-in Glen's flow law (TODO).
Newton linearisation of the viscosity can be used using the keywords Nonlinear System Newton After Tolerance and/or Nonlinear System Newton After Iterations. It is automatically reset to False at the beginning of a new iteration.

The Mean Density and Mean Viscosity, if not uniform along the vertical direction, can be computed using the GetMeanValueSolver routine or the StructuredProjectToPlan solver (preferred solution).

Contrary to the NS solver, the gravity must be orientated along the z-axis and is taken from the value of Flow BodyForce 2 for a SSA-1D problem or Flow BodyForce 3 for a SSA-2D problem.

A Neumann condition on the lateral boundaries can be applied with the keyword Calving front = Logical True in the Boundary condition section. The condition is : *0.5 * g * (rho_ice * h^2 - rho_water * h_im^2)* where
- *g* is the absolute value of the gravity taken from Flow BodyForce i
- *rho_ice* is the ice Mean Density
- *rho_water* is water density taken from the constants section (or default=1.03225e-18)
- *h* is the front thickness computed as Zs-Zb
- *h_im* is the thickness below sea level computed as Sea Level - Zb, where Sea Level is taken from the constants section (or default=0.0).
Note that in the absence of explicit boundary condition (no dirichlet condition or Calving front = Logical True not found) the natural boundary condition is force equilibrium *(rho_ice * h^2 = rho_water * h_im^2)*.

The SSA velocities and pressure can be used, for example, as initial conditions for the Stokes Solver.

When the SSA solution is computed on a boundary of a mesh of dimension larger than the SSA problem (e.g. a 3D mesh for a SSA-2D problem), the SSA solution computed on the boundary can be:
- exported on the whole mesh using the StructuredProjectToPlane solver (preferred solution) or the SSASolver routine
- used as a Dirichlet condition for the SIA velocity (see the [SIA Solver](./SIA.md)).

### Basal friction

#### friction laws
Currently, there are 4 friction laws implemented in the SSA solver:

- a linear friction law
*tau_b = beta . u*
- a Weertman type friction law
*tau_b = beta.{u_b}^{m - 1} . u*
- a Coulomb type friction law
*tau_b = 1/{A_s}^{1/n} {[{ 1/ {(1 + alpha . chi^q)} }]}^{1/n} . {u_b}^{1/n-1}. u*
where *alpha = {(q - 1)^{q-1}}/{q^q}* and *chi = {u_b}/{C^n N^n A_s}*
- a regularised coulomb friction law  
*tau_b = beta (u_b/(u_b+u_0))^m*

The latter three are non-linear and a Newton linearisation can be used. 
When *u_b = (u^2+v^2)^{1/2}< u_min*, *u_b* in the previous equations is replaced by *u_min*.

The friction law is chosen using the keyword SSA Friction Law, which takes the value Linear, Weertman, coulomb, regularised Coulomb. The other keywords are:  
- a linear friction law
  - SSA Friction Parameter → *beta*
- a Weertman type friction law
  - SSA Friction Parameter → *beta*
  - SSA Friction Exponent → *m*
  - SSA Friction Linear Velocity → *u_lin*
- a Regularised Coulomb friction law  
  - SSA Friction Parameter → *beta*  
  - SSA Friction Exponent → *m*
  - SSA Friction Linear Velocity → *u_lin*
  - SSA Friction Threshold Velocity → *u_0*
- a Coulomb type friction law
  - SSA Friction Parameter → *beta= {A_s}^{-m}*
  - SSA Friction Exponent → *m = 1/n*
  - SSA Friction Linear Velocity → *u_lin*
  - SSA Friction Post-Peak → *q >= 1*
  - SSA Friction Maximum Value → *C ~ max bed slope*
  - Effective Pressure (variable) → *N*
  - SSA Min Effective Pressure → *N_{min}*, such that *N >= N_{min}*

#### Sub-Element grounding line parametrisation
The flotation condition can be tested directly at the integration points. The friction parameter is then set to 0 if ice is floating at the integration point.
This scheme is activated with the Solver Keyword *Sub-Element GL parameterization = logical True*. The number of integration point for partially floating elements, i.e. where the scheme is used is set with the keyword
*GL integration points number = N*

:warning: if this scheme is not used it is to the user responsibility to parameterised *beta* a a function of the grounded mask and set it to 0 for floating nodes.

#### Post-processing variables

- if a variable named "strbasemag" is found computes the element-averaged friction
- if a variable names "Ceff" is found computes the nodal effictive friction coefficient

## SIF contents
Solver section:

```
Solver 1
  Equation = "SSA"
  Procedure = File "ElmerIceSolvers" "SSABasalSolver"
  Variable = String "SSAVelocity"
  Variable DOFs = 2   ! 1 in SSA 1-D or 2 in SSA-2D

  !# Set true to set friction to 0 when ice is floating
  !   Require: GroundedMask and Bedrork variables
  Sub-Element GL parameterization = logical True
  !# You can set the number of IP points to evaluated the friction in the
  !   first floating element
  GL integration points number = Integer 20

  Linear System Solver = Direct
  Linear System Direct Method = umfpack

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-08
  Nonlinear System Newton After Iterations = 5
  Nonlinear System Newton After Tolerance = 1.0e-05

  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = Real 1.0e-3
End
Material Properties:

Material 1

! Material properties
  Viscosity Exponent = Real $1.0/n
  Critical Shear Rate = Real 1.0e-10

  SSA Mean Viscosity = Real $eta
  SSA Mean Density = Real $rhoi

! Needed for Linear, Weertman , Coulomb and regularised coulomb
  ! Which law are we using (linear, weertman , coulomb or regularised coulomb)
  SSA Friction Law = String "Coulomb"
  ! beta parameter (beta = 1/As^m)
  SSA Friction Parameter = Variable coordinate 1 , Coordinate 2
     Real  MATC "1.0e-3*(1.0 + sin(2.0*pi* tx(0) / L)*sin(2.0*pi* tx(1) / L))

! Needed for Weertman , Coulomb and regularised coulomb
  ! Exponent m 
  SSA Friction Exponent = Real $1.0/n
  
  ! Min velocity for linearisation where ub=0
  SSA Friction Linear Velocity = Real 0.0001

! Needed for Coulomb only
  ! post peak exponent in the Coulomb law (q, in Gagliardini et al., 2007)
  SSA Friction Post-Peak = Real 1.0
  ! Iken's bound  tau_b/N < C (see Gagliardini et al., 2007)
  SSA Friction Maximum Value = Real 0.5

! Needed for Regularised Coulomb
  SSA Friction Threshold Velocity = Real 300.0

End
Body Forces:

Body Force 1
  Flow BodyForce 1 = Real 0.0
  Flow BodyForce 2 = Real 0.0
  Flow BodyForce 3 = Real $gravity
End
Constants:

Constants
! Used for Neumann condition
  Water Density = Real ....
  Sea Level = Real ...
End
Boundary Conditions:

Boundary Condition 1 
! Dirichlet condition
  SSAVelocity 1 = Real ...
  SSAVelocity 2 = Real ...
End
Boundary Condition 1 
! Neumann Condition
  Calving Front = Logical True
End
```
For the “GetMeanValueSolver” routine, the required keywords in the SIF file for this solver are:

```
Solver 1
  Equation = "SSA-IntValue"
  Procedure = File "ElmerIceSolvers" "GetMeanValueSolver"
  Variable = -nooutput String "Integrated variable"
  Variable DOFs = 1

  Exported Variable 1 = String "Mean Viscosity"
  Exported Variable 1 DOFs = 1
  Exported Variable 2 = String "Mean Density"
  Exported Variable 2 DOFs = 1

  Linear System Solver = Direct
  Linear System Direct Method = umfpack

  Steady State Convergence Tolerance = Real 1.0e-3
End

!!! Upper free surface
Boundary Condition 1

  Depth = Real 0.0
  Mean Viscosity = Real 0.0
  Mean Density = real 0.0
End
```
For the “SSASolver” routine, the required keywords in the SIF file for this solver are:

```
Solver 4
  Equation = "SSA Velocity"
  Procedure = File "ElmerIceSolvers" "SSASolver"
  Variable = -nooutput String "varSSA"
  Variable DOFs = 1

  Exported Variable 1 = String "SSAFlow"
  Exported Variable 1 DOFs = 4  ! 3 in 2D, 4 in 3D

  Linear System Solver = Direct
  Linear System Direct Method = umfpack

  Steady State Convergence Tolerance = Real 1.0e-3
End

!!! bedrock
Boundary Condition 1
  
  SSAFlow 1 = Equals SSAVelocity 1
  SSAFlow 2 = Equals SSAVelocity 2
  SSAFlow 3 = Real 0.0e0
End

!!! Upper free surface
Boundary Condition 2
 
  Depth = Real 0.0
 
  SSAFlow 4 = Real 0.0     ! p=0 at the surface
End
```
If one wants to solve the SSA + SIA, the sif will read:

```
Solver 4
  Equation = "SIA Velocity"
  Procedure = File "SIASolver" "SIASolver"
  Variable = -nooutput String "varSIA"
  Variable DOFs = 1

  Exported Variable 1 = String "SIAFlow"
  Exported Variable 1 DOFs = 4  ! 3 in 2D, 4 in 3D

  Linear System Solver = Direct
  Linear System Direct Method = umfpack

  Steady State Convergence Tolerance = Real 1.0e-3
End

!!! bedrock
Boundary Condition 1
  
 ...
  SIAFlow 1 = Equals SSAVelocity 1
  SIAFlow 2 = Equals SSAVelocity 2
  SIAFlow 3 = Real 0.0e0
  ...
End

!!! Upper free surface
Boundary Condition 2
 
 ...  
  SIAFlow 4 = Real 0.0     ! p=0 at the bottom
End
```

## Examples
For examples look in your elmer source distribution under [ELMER_TRUNK]/elmerice/Tests/SSA and under [ELMER_TRUNK]/elmerice/examples/Test_SSA.
