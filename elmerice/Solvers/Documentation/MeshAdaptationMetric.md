# Mesh Adaptation (metric)
## General Information
- **Solver Fortran File:** MMG2D_MetricAniso.F90
- **Solver Name:** ElmerIce_MeshAdapt2D(MMG2D_MetricAniso)
- **Required Output Variable(s):**
  - (1) Metric (dofs = 3)
  - (2) hessian (dofs = 3)
- **Required Input Variable(s):**
  - (1) Nodal gradient (dofs = 2)
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This solver is used for the mesh adaptation ([Mesh Adaptation](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:meshadaptation)) to compute the anisotropic metric M.

The metric M , used to define the element size, derives from a geometric error estimate based on an upper bound for the interpolation error of a continuous field to piecewise linear elements (Frey and Alauzet, 2005).

For a variable v, M depends on the eigenvalues lambda_i and eigenvector matrix R of the hessian matrix of v, H (i.e. small elements are required where the curvature is the highest):

*M=R.Lambda.R^{-1} with Lambda=(matrix{2}{2}{lambda_1 0  0 lambda_2}) and \lambda_i=min ( max ( {{c|\lambda_i|}/{epsilon_v}},{{1}/{l^2_{max}}} ), {{1}/{l^2_min}} ),    (1)*

where

- *c* is a geometric constant equal to 2/9 in 2D
- *l_min (resp. l_max)* is a prescribed minimal (resp. maximal) edge size
- *epsilon_v* is the prescribed maximum error
First this solver computes the hessian matrix *H*. As computing second derivatives in linear elements is not straightforward this is done by solving the diffusive equation *H_ij+ K ∇(H_ij) = 1/2 (dq_i/dx_j+dq_j/dx_i)*, where *K=k A* is a diffusivity proportionnal to the local element size (*A*) and *g_i={{\partial v}/{\partial x_i}}* are the nodal gradients of the variable *v* (This can be computed using using the [Compute2DNodalGradient](./NodalGradient.md) Solver).

Finally, the metric *M* is then computed from Eq. (1)

## SIF contents
```
Solver 5
   Equation = "Metric2"
   Variable = -nooutput dumy

   Procedure = "ElmerIce_MeshAdapt2D" "MMG2D_MetricAniso"

   Metric Variable Name = String "M2"

   Hessian Variable Name = String "ddx2"
   Gradient Name = String "Gradient2"
   Diffusivity = Real 0.5 !! the diffusivity k; the total diffusivity is kA

   Linear System Solver = Direct
   Linear System Direct Method = umfpack

  Exported Variable 1 = -dofs 3 "M2"
  Exported Variable 2 = -dofs 3 "ddx2"
End


Body Force 1
!! Parameters in Eq. 1
  M2 Hmin = Real 1.0e-3
  M2 Hmax = Real 1.0
  M2 err =  Real 0.0033
End
```

## Example
Examples for anisotropic mesh adaptation can be found under [ELMER_TRUNK]/elmerice/Tests/MMG2D_Aniso1 and [ELMER_TRUNK]/elmerice/Tests/MMG2D_Aniso2, where the mesh size is adapted using 1 or 2 variables (i.e. combining metric information), respectively.
