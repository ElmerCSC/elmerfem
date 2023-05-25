## AdjointThickness : Direct Solver {#thickness_direct_solver}

**Module name**: AdjointThickness_ThicknessSolver.F90      
**Module subroutines**: AdjointThickness_ThicknessSolver     
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)     
**Document authors**: Fabien Gillet-Chaulet   
**Document edited**: 10/12/2020   

**Solver Variable:**

 - H : thickness


### Introduction

This solver solves the steady-state Thickness evolution equation:
$$div(uH) = Ms + Mb$$
where:    
- $u$ is the mean horizontal velocity    
- $Ms$ and $Mb$ are the surface and bottom mass balance (>0 for accumulation)    

The convection velocity can be provided as a variable using the keyword *Flow Solution Name* in the *Solver Section*, or **alternatively** the velocity componenents can be read in the *body force section.* 

In general $Ms$ and $Mb$ will be the *apparent mass basalance*, i.e. they will include a correction with the observed thickness rate of change.

This solver has been derived from the legacy Elmer/Ice Thickness solver. It has been simplified and separated from the legacy solver to derive the adjoint code ([AdjointThickness : Gradient Solver](#thickness_gradient_solver))
and make sure that they are consistent.

In particular, this solver allows to use only the **stabilised method** (i.e. a StreamLine Upwind Petrov-Galerkin (SUPG) method) and **does not include the limiters**.

It can only be used in a **cartesian reference frame** and in **steady-state**.

Note that in the absence of limiters this equation is linear and does not require non-linear iterations.

### Keywords


#### Solver Section:

```
  Solver *id*
   Equation = "Thickness"

  !# The solver variable 
   Variable = "H"

   Procedure = "ElmerIceSolvers" "AdjointThickness_ThicknessSolver"

 !# Keywords related to linear system
  Linear System Solver = ...

 !#  the convection velocity $u$
  Flow Solution Name = String  ...

 !# if the convection velocity is read from the body force, 
 !#  provide the diemnsion
 Convection Dimension = Integer
End
```

#### Body Forces:

```
!# Ms and Mb
Body Force *id*
  Top Surface Accumulation = Real ...
  Bottom Surface Accumulation = Real ...

  !# If <Flow Solution Name> is not provided in the solver section
  !#  the convection velocity can be read here:
  ! x-componenent
  Convection Velocity 1 = Real ...
  ! y-componenent if Convection Dimension = 2
  Convection Velocity 1 = Real ...
End
````

#### Boundary Conditions:

Requires Dirichlet boundary conditions at inflow boundaries.
```
Boundary Condition 1 
! Dirichlet condition
  H = Real ...
End
```

### Tests and Examples

- See examples for the [Mass Conservation methods](../../examples/Inverse_Methods)
