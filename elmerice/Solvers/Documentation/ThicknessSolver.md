#Solver for the Thickness evolution equation
##General Information
- **Solver Fortran File:** ThicknessSolver.f90
- **Solver Name:** ThicknessSolver
- **Required Output Variable(s):** H
- **Required Input Variable(s):** H residual
- **Optional Output Variable(s):** dhdt
- **Optional Input Variable(s):** FlowSolution

##General Description
Solve the Thickness evolution equation:
*dh/dt + div(uH) = Ms + Mb*
where:
- *u* is the mean horizontal velocity
- *Ms* and *Mb* are the surface and bottom mass balance (>0 for accumulation)
This solver is based on the FreeSurfaceSolver and use a SUPG stabilsation scheme by default (residual free bubble stabilization can be use instead).
As for the Free surface solver Min and Max limiters can be used.
As for the Free surface solver only a Dirichlet boundary condition can be imposed.
This solver can be used on a mesh of the same dimension as the problem (e.g. solve on the bottom or top boundary of a 3D mesh to solve the 2D thickness field) or on a mesh of lower dimension (e.g. can be use in a 2D plane view mesh with the [SSA solver](./SSA.md) for example).

When working on a mesh of the same dimension as the problem it can be usefull to have an extruded mesh along the vertical direction and to use the StructuredProjectToPlane and StructuredMeshMapper solvers to compute the mean horizontal velocity (from the Stokes solution), export the value of H computed on one boundary in the whole mesh and update the mesh (see examples).

##SIF contents
Solver section:

```
Solver 1
   Equation = "Thickness"
   Variable = -dofs 1 "H"

   Exported Variable 1 = -dofs 1 "H Residual"

!! To compute dh/dt
   Exported Variable 2 = -dofs 1 "dHdt"
   Compute dHdT = Logical True

  Procedure = "ElmerIceSolvers" "ThicknessSolver"
   Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

   Linear System Solver = Iterative
   Linear System Max Iterations = 1500
   Linear System Iterative Method = BiCGStab
   Linear System Preconditioning = ILU0
   Linear System Convergence Tolerance = Real 1.0e-12
   Linear System Abort Not Converged = False
   Linear System Residual Output = 1500

! equation is linear if no min/max
   Nonlinear System Max Iterations = 50
   Nonlinear System Convergence Tolerance  = 1.0e-6
   Nonlinear System Relaxation Factor = 1.00

! stabilisation method: [stabilized\bubbles]
  Stabilization Method = stabilized
  
!! to apply Min/Max limiters
  Apply Dirichlet = Logical True

!! to use horizontal ALE formulation
   ALE Formulation = Logical True

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
!! Limiters
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

##Examples
For examples look in your elmer source distribution under [ELMER_TRUNK]/elmerice/Tests/ThicknessSolver or [ELMER_TRUNK]/elmerice/Tests/SSA_IceSheet or [ELMER_TRUNK]/elmerice/examples/Test_SSA.
