# Solver m1qn3
## General Information
- **Solver Fortran File:** m1qn3.f
- **Solver Name:** Optimize_m1qn3

## General Description
[m1qn3.pdf](./m1qn3.pdf)

## SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 6
  Equation = "Optimize_m1qn3"

  procedure = "ElmerIceSolvers" "Optimize_m1qn3Parallel"

  Cost Variable Name = String "CostValue"
  Optimized Variable Name = String "Beta"
  Gradient Variable Name = String "DJDBeta"
  gradient Norm File = String "GradientNormAdjoint_$name$.dat"
  Mesh Independent = Logical True

! M1QN3 Parameters
  M1QN3 dxmin = Real 1.0e-10
  M1QN3 epsg = Real  1.e-6
  M1QN3 niter = Integer 200
  M1QN3 nsim = Integer 200
  M1QN3 impres = Integer 5
  M1QN3 DIS Mode = Logical False
  M1QN3 df1 = Real 0.5
  M1QN3 normtype = String "dfn"
  M1QN3 OutputFile = File  "M1QN3_$name$.out"
  M1QN3 ndz = Integer 20
```

## Examples
A 3D examples can be found in [ELMER_TRUNK]/elmerice/Tests/InvMeth*.

## Note about Continuous/Discrete Gradients
It is recommended to use the 'Mesh Independent = Logical True' keyword, as this accounts for the variable 'footprint' (i.e. element area) of different nodes on the base of the mesh. This issue becomes particularly noticeable when there is a large variation in mesh resolution across the domain. In such cases, without 'Mesh Independent = Logical True', M1QN3 will tend to alter Beta values only in regions where elements are large (as these make the largest contribution to the discrete cost gradient).

## Solver's references to cite
Since its version 3.2 (December 2008), M1QN3 is distributed under the GNU General Public License and listed in the Free Software Directory. You are also kindly asked to mention the following paper in your publications for which M1QN3 has been used:

J.Ch. Gilbert and C. Lemaréchal (1989). Some numerical experiments with variable-storage quasi-Newton algorithms. Mathematical Programming 45, 407-435.
