## AdjointSSA : Direct Solver {#ssa_direct_solver}

**Module name**: AdjointSSA_SSASolver.F90  
**Module subroutines**: AdjointSSA_SSASolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)    
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 28/04/2020

**Solver Variable:**

 - SSAVelocity : horizontal velocity  

**Required input variables:**

  - Zb: bottom surface elevation
  - Zs: top surface elevation
  - *[Optional - see below]* GroundedMask : mask for grounded (>=0) or floating (<0) ice
  - *[Optional - see below]* bedrock : bed elevation


### Introduction

This solver solves the classical SSA equation for the horizontal velocity.

It is very similar to the legacy Elmer/Ice [SSABasalSolver](http://elmerfem.org/elmerice/wiki/doku.php?id=solvers:ssa).
However it is separated to derive the adjoint code ([AdjointSSA : Gradient Solver](#ssa_gradient_solver))
and make sure that they are consistent.


### Keywords


#### Solver Section:

```
  Solver *id*
  Equation = "SSA"

  !# dimension 1 or 2 for flowline or 2D
   Variable = -dofs 2 "SSAVelocity"

   Procedure = "ElmerIceSolvers" "AdjointSSA_SSASolver"

  !# Set true to set friction to 0 when ice is floating
  !   Require: GroundedMask and Bedrork variables
  Sub-Element GL parameterization = logical True
  !# You can set the number of IP points to evaluated the friction in the
  !   first floating elemnt
  GL integration points number = Integer 20

  !# Keywords related to linear system
  Linear System Solver = ...

  !# Keywords related to non-linear iterations 
  !  for the accuracy of the gradient computation
  !   Newton has to be used if the problem is non linear
   Nonlinear System Max Iterations = Integer ...
   Nonlinear System Convergence Tolerance  = Real ....
   Nonlinear System Newton After Iterations = Integer ...
   Nonlinear System Newton After Tolerance = Real ...
   Nonlinear System Relaxation Factor = 1.00


  End

```

#### Material Properties:
```
  Material *id*

  Viscosity Exponent = Real ....
  Critical Shear Rate = Real ....

  SSA Mean Viscosity = Real ...
  SSA Mean Density = Real ...

 !# Friction law (linear or weertamn)
  SSA Friction Law = String ...
 !# Friction coefficient
 !  it is set to 0 automatically for floating ice if 
 !    Sub-Element GL parameterization = logical True [see above]
  SSA Friction Parameter = Real ...

 !# Keywords related to Weertman
 !# SSA Friction Exponent m 
  SSA Friction Exponent = Real
 ! Min velocity for linearisation where ub=0
  SSA Friction Linear Velocity = Real ....

```

#### Constants: 
```
Constants
! Used for Neumann condition
  Water Density = Real ....
  Sea Level = Real ...
End
```

#### Body Forces:
```
!# gravity for BodyForce 2 (1D-case) or BodyForce 3 (2D plane view case)
Body Force *id*
  Flow BodyForce 1 = Real 0.0
  Flow BodyForce 2 = Real 0.0
  Flow BodyForce 3 = Real $gravity
End
````

#### Boundary Conditions:
```
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

### Tests and Examples

- See examples for the [SSA inverse methods](../../examples/Inverse_Methods)
