# Mesh Adaptation (metric intersection)
## General Information
- **Solver Fortran File:** MMG2D_MetricIntersect.F90
- **Solver Name:** ElmerIce_MeshAdapt2D(MMG2D_MetricIntersect)
- **Required Output Variable(s):**
  - (1) Metric 1_2
- **Required Input Variable(s):**
  - (1) Metric 1
  - (2) Metric 2
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This is a pseudo solver (i.e. it is not solving an equation).

This solver is used for the mesh adaptation ([Mesh Adaptation](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:meshadaptation)) to perform the intersection of two anisotropic metrics *M_1* and *M_2*.

The intersection *M_{1∩2}=M_1 ∩ M_2* of two metrics *M_1* and *M_2* is given by (Alauzet et al., 2007):

*M_{1∩2}=^{T}P^{-1} (matrix{2}{2}{{max(mu^1_{1},mu^2_{1})} 0  0 {max(mu^1_{2},mu^2_{2})}})P^{-1}*

with *P* the matrix where the columns are the normalised eigenvectors *(e_i)_{i=1,2}*, of *N=M^{-1}_{1}M_2* and *μ^j_{i}=^{T}e_{i}M_{j}e_i*.

*M_1* and *M_2* can be computed using the [MMG2D_MetricAniso](./MeshAdaptationMMG2D.md) Solver.

A variable containing the metric *M_i* must have 3 dofs *(M_{11},M_{22},M_{12})*.

## SIF contents
```
Solver 6
   Equation = "Metric"
   Variable = -nooutput dumy

   Procedure = "ElmerIce_MeshAdapt2D" "MMG2D_MetricIntersect"

  Optimize Bandwidth = False


   Metric Variable Name = String "M1M2"
   Metric 1 Variable Name = String "M2"
   Metric 2 Variable Name = String "M1"

   Exported Variable 1 = -dofs 3 "M1M2"
End
```

## Example
An example for anisotropic mesh adaptation using 2 variables can be found under [ELMER_TRUNK]/elmerice/Tests/MMG2D_Aniso2.
