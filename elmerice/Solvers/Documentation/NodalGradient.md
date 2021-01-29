# Nodal Gradient
## General Information
- **Solver Fortran File:** Compute2DNodalGradient.F90
- **Solver Name:** ElmerIce_MeshAdapt2D(Compute2DNodalGradient)
- **Required Output Variable(s):**
  - (1) g (dofs=2)
- **Required Input Variable(s):**
  - (1) v
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This is a pseudo solver (i.e. it is not solving an equation). This solver compute the nodal 2D gradient vector *g_i={{\partial v}/{\partial x_i}}* of a variable *v*.

This is used for example by the mesh adaptation procedure to compute the hessian of *v* ([Mesh Adaptation](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:meshadaptation)).

By default (FE consistent average = Logical True), this is done using a L²-projection on the FE mesh; If FE consistent average = False, at a given node the derivative is simply the average of the derivatives evaluated at the node in each elements sharing the node.

## SIF contents
```
Solver 2
  Equation = "Nodal Gradient"
  Variable = -dofs 2 "g"
  Procedure = "ElmerIce_MeshAdapt2D" "Compute2DNodalGradient"

  Optimize Bandwidth = False

  Variable Name = string "v"
  FE consistent average = Logical True
End
```

## Example
Examples for anisotropic mesh adaptation can be found under [ELMER_TRUNK]/elmerice/Tests/MMG2D_Aniso1 and [ELMER_TRUNK]/elmerice/Tests/MMG2D_Aniso2, where the mesh size is adapted using 1 or 2 variables (i.e. combining metric informations), respectively.
