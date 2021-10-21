# Solver for the Thickness evolution equation
## General Information
- **Solver Fortran File:** ThicknessSolver.f90  
- **Solver Name:** ThicknessSolver  
- **Required Output Variable(s):** H  
- **Required Input Variable(s):** H residual  
- **Optional Output Variable(s):** dhdt  
- **Optional Input Variable(s):** FlowSolution  

## General Description

### Description

Solve the Thickness evolution equation:  
$\dfrac{dH}{dt} + div(uH) = M_s + M_b$  
where:  
- $u$ is the mean horizontal velocity  
- $M_s$ and $M_b$ are the (vertical) surface and bottom mass balances (>0 for accumulation)

This solver is based on the FreeSurfaceSolver and use a SUPG stabilisation scheme by default (residual free bubble stabilization can be use instead).  By default the SUPG stabilisation parameter is computed as  
$\tau_{SUPG}^{stead}= h_k / 2 ||u||$ with $h_k$ the element diameter and ||u|| the norm of the convection velocity $u$.  
Several authors have noted that $\tau_{SUPG}$ should be limited by $\Delta t / 2$ in transient (see e.g. Akin and Tezduyar, 2004). Here if the flag  *Transient Stabilisation* is set to **TRUE** the stabilisation parameter is computed as:  
$\tau_{SUPG}=( 1/\tau_1^r + \alpha 1/\tau_2^r)^{(-1/r)}$   with $r=2$ and, similarly to $\tau_{SUPG}^{stead}$, $\tau_1$ estimate the advection time through the element as   $\tau_1=(\sum_{a=1}^{n_a} |u. \nabla \Phi_a)^{-1}$ with $\Phi_a$ the basis function attached to node $a$  and $\tau_2=\Delta t / 2$.   

Here $\alpha=1.0e-3$ and can be set with the flag *Tau2 factor = Real ...*. Setting $\alpha=0.0$ leads to the steady state limit, i.e. $\tau_{SUPG}=\tau_1$, and an upper limit $\alpha=1.0$ leads to the transient limit $\tau_{SUPG}=\Delta t / 2$ for large advetion times ($\tau_1$).


As for the Free surface solver Min and Max limiters can be used. 

As for the Free surface solver only  Dirichlet boundary conditions can be imposed.  

This solver can be used on a mesh of the same dimension as the problem (e.g. solve on the bottom or top boundary of a 3D mesh to solve the 2D thickness field) or on a mesh of lower dimension (e.g. can be use in a 2D plane view mesh with the [SSA solver](./SSA.md) for example).

When working on a mesh of the same dimension as the problem it can be useful to have an extruded mesh along the vertical direction and to use the StructuredProjectToPlane and StructuredMeshMapper solvers to compute the mean horizontal velocity (from the Stokes solution), export the value of H computed on one boundary in the whole mesh and update the mesh (see examples).

### Remarks  
- the equation is linear and becomes non-linear (i.e. requires non-linear iterations only if limiters are used).

### References

Akin, J.E., Tezduyar, T.E., 2004. Calculation of the advective limit of the SUPG stabilization parameter for linear and higher-order elements. Computer Methods in Applied Mechanics and Engineering 193, 1909–1922. https://doi.org/10.1016/j.cma.2003.12.050


### Possible evolutions  
- Limiters are now a generic Elmer functionality so it might be better to use this functionality instead of the internal limiters (cf [ElmerSolver Manual](http://www.nic.funet.fi/pub/sci/physics/elmer/doc/ElmerSolverManual.pdf))
- If requested the *dhdt* is computed as the first order time differential of H; This is also now a generic Elmer feature that can be activated with the keyword *Calculate Velocity* (see section **13.4 Exported and derived variables** in the [ElmerSolver Manual](http://www.nic.funet.fi/pub/sci/physics/elmer/doc/ElmerSolverManual.pdf))

## SIF contents
Solver section:

```
Solver 1
   Equation = "Thickness"
   Procedure = "ElmerIceSolvers" "ThicknessSolver"
   Variable = -dofs 1 "H"

!! Numerical settings, e.g.:
 !!! Linear system resolution
   Linear System Solver = Iterative
   Linear System Max Iterations = 1500
   Linear System Iterative Method = BiCGStab
   Linear System Preconditioning = ILU0
   Linear System Abort Not Converged = False
   Linear System Residual Output = 1500

  !! Mandatory if Apply Dirichlet = True 
   Linear System Convergence Tolerance = Real 1.0e-12

 !!! Non linear system resolution
 !!! Rq. equation is noon linear only if limiters are used
   Nonlinear System Max Iterations = 50
   Nonlinear System Convergence Tolerance  = 1.0e-6
   Nonlinear System Relaxation Factor = 1.00

!!! stabilisation method: [stabilized/bubbles] 
  Stabilization Method = stabilized ![Default]

!!! Activate Transient Stabilisation scheme  
  Transient Stabilisation = Logical True ![Default: False]
  Tau2 factor = Real ... ![Default: 1.0e-3]

!! To compute dh/dt
   Compute dHdT = Logical True ![Default: False]
   Exported Variable 2 = -dofs 1 "dHdt"
  
!! to apply Min/Max limiters
  Apply Dirichlet = Logical True ![Default: False]
 !! Mandatory if using internal limiters
  Exported Variable 1 = -dofs 1 "H Residual"

!! to use horizontal ALE formulation
   ALE Formulation = Logical True ![Default: False]

!! To get the mean horizontal velocity

!!  either give the name of the variable
     Flow Solution Name = String "SSAVelocity"
     
!!!!! or give the dimension of the problem using:
! Convection Dimension = Integer
End
```
Material Properties:

```
Material 1
!! Internal Limiters
  Min H = Real ....
  Max H = Real ....

End
```
Body Forces:

```
Body Force 1
!! Mass balance
  Top Surface Accumulation = Real ....
  Bottom Surface Accumulation = Real ....
  
  
!! if the convection velocity is not directly given by a variable
!! Then give //Convection Dimension = Integer// in the solver section 
!! and the Mean velocity here:

  Convection Velocity 1 = Variable int Velocity 1, thickness
     REAL MATC "tx(0)/tx(1)"
  Convection Velocity 2 = Variable int Velocity 2, thickness
     REAL MATC "tx(0)/tx(1)"
 
 
End
```
Boundary Conditions:

```
Boundary Condition 1 
! Dirichlet condition only
  H = Real ...
End
```

## Examples
For examples look in your elmer source distribution under [ELMER_TRUNK]/elmerice/Tests/ThicknessSolver or [ELMER_TRUNK]/elmerice/Tests/SSA_IceSheet or [ELMER_TRUNK]/elmerice/examples/Test_SSA.
