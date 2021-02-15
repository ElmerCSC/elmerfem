# Solver FlowDepth Solver
## General Information
- **Solver Fortran File:** Flowdepth.f90
- **Solver Name:** FlowdepthSolver
- **Required Output Variable(s):** Depth or Height (user defined)
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** FreeSurf, FreeSurfGrad1, FreeSurfGrad2 (FreeSurf name is user defined)
- **Optional Input Variable(s):** None

## General Description
This solver computes the vertical distance from a surface (either free surface or bedrock in Glaciology). In structured FD-grids this is a trivial procedure - not so in unstructured FEM meshes as used by Elmer. A degenerated Poisson equation has to be solved in order to get nodal values for either the flow depth below the free surface or the flow height above the bedrock. It can also compute the horizontal derivative of the upper/lower surface if Depth/height is calculated.

An optional feature provided by this solver is to compute the gradient components of the surface from where the distance is evaluated.

Alternatively, if you are using a structured mesh, you might consider replacing the FlowDepth solver by the StructuredProjectToPlane solver, as explained [here](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:structuredmesh).

## SIF contents
The required keywords in the SIF file for this solver are:

```
! This solves the depth underneath the free surface
!--------------------------------------------------
Solver 1
  Equation = "Flowdepth"
   Exec Solver = "Before Timestep"
   Procedure = File "ElmerIceSolvers" "FlowDepthSolver"
   Variable = String "Depth"
   Variable DOFs = 1
   Linear System Solver = "Direct"
   Linear System Direct Method = "UMFPACK"
   ! this sets the direction
   ! -1 is negative z-direction (upside down)
   ! +1 is positive (downside up)
   Gradient = Real -1.0E00
  ! switch that to True, if you want to have 
  ! free surface gradients to be computed
  !------------------------------------
  Calc Free Surface = Logical True
  ! the name for the exported (if not existing) added variable
  ! the gradients will be stored in variables with the base
  ! name given and "Grad1" and (in 3 dimensions) "Grad2" added,
  ! so in our case "FreeSurfGrad1" and "FreeSurfGrad2"
  ! again, if those variables did not exist, they will be
  ! automatically created
  !-----------------------------------------------------------
  Freesurf Name = String "FreeSurf"
End

! This solves the height above the bedrock
!------------------------------------------
Solver 2
  Equation = "Flowheight" ! mind different name
  Exec Solver = "Before Timestep"
  Procedure = File "Flowdepth2" "FlowDepthSolver" ! make a copy of the original solver, to get a separate address space
  Variable = String "Height" ! mind different name for variable
  Variable DOFs = 1
  Linear System Solver = "Direct"
  Linear System Direct Method = "UMFPACK"
  
  Gradient = Real 1.0E00 ! this time positive
  Calc Free Surface = Logical False
End

! boundary at the free surface
!-----------------------------
Boundary Condition 1
  ...
  Depth = Real 0.0
End

! boundary at bedrock
!-----------------------------
Boundary Condition 2
  ...
  Height = Real 0.0
End

! internal boundary with
! no Dirichlet and neither
! the default Neumann condition
!-----------------------------
Boundary Condition 3
  ...
  Flowdepth Skip = Logical True
End
```
## Examples
Tests using the FlowDepth and FlowHeight solvers can be found in [ELMER_TRUNK]/elmerice/Tests/FlowDepth and [ELMER_TRUNK]/elmerice/Tests/FlowHeight .


