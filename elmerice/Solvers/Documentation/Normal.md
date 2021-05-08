# Solver ComputeNormal
## General Information
- **Solver Fortran File:** ComputeNormal.f90
- **Solver Name:** ComputeNormal
- **Required Output Variable(s):** default is Normal Vector
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
Solver for computing nodal normal vectors on boundaries based on an averaged value of the adjacent elements. Bear in mind that in parallel computations this requires halo-elements, such that cross-partition information is available.

## SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 1
   Exec Solver = "Before Simulation"
   Equation = "Normal vector"
   Variable = "Normal Vector"   
   ! in 3dimensional simulations we have 3 entries
   Variable DOFs = 3 
   !NB: does not need to actually solve a matrix
   !    hence no BW optimization needed
   Optimize Bandwidth = Logical False 
   Procedure = "ElmerIceSolvers" "ComputeNormalSolver"
   ! if set to True, all boundary normals would be computed by default
   ComputeAll = Logical False
End

! on this boundary, we want the normals to be computed
Boundary Condition 1
    ComputeNormal = Logical True
End

! on this boundary, we want to skip computation of the normals
!       (default, but overrulled by ComputeAll = Logical True)
Boundary Condition 2
    ComputeNormal = Logical False
End
```

## Examples
A 2D example can be found in [ELMER_TRUNK]/elmerice/Tests/Contact.
