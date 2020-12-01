# Mesh Adaptation (MMG2D)
## General Information
- **Solver Fortran File:** MMG2DSolver.F90
- **Solver Name:** ElmerIce_MeshAdapt2D(MMG2DSolver)
- **Required Output Variable(s):**
  - (1) dumy
- **Required Input Variable(s):**
  - (1) Metric
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
This is a pseudo solver (i.e. it is not solving an equation). This solver is used for the mesh adaptation ([Mesh Adaptation](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:meshadaptation)).

This solver call the Mmg library (http://www.mmgtools.org/) to perform the mesh adaptation. This requires you to separately install the Mmg code (https://github.com/MmgTools/mmg). This solver is included only if the cmake arguments “MMG_INCLUDE_DIR” and “MMG_LIBRARY” are provided when compiling elmerice.

It will perform isotropic or anisotropic 2D mesh adaptation depending on the dofs of the Input Variable:

- If dofs=1, it must contain the required element size at the current node location
- If dofs=3, it must contain the components (M11,M22,M12) of the anisotropic metric (see [Mesh Adaptation](http://elmerfem.org/elmerice/wiki/doku.php?id=mesh:meshadaptation))
Our implementation is actually restricted to the adaptation of plane-view 2D meshes comprised of linear 3-nodes triangular elements.

Tested with Mmg master branch commit (6acfa9e7b20e41134d56af10eba1bb8fd1283f8f).

## SIF contents
```
Solver 1
  Exec Solver = after timestep

  Equation = "MMG"
  Variable = -nooutput dumy
    Procedure = "ElmerIce_MeshAdapt2D" "MMG2DSolver"

!! Name of the adapted mesh (will be used to save the mesh on disk)
  Output file name = "square_aniso"

!! Name of the variable that contain the metric; 
  !! Anisotropic 2D mesh adaptation if M is of size 3 (M11,M22,M12)
  !! Isotropic 2D mesh adaptation if M is of size 1 
  Metric Variable Name = String "M" 
  

!! Mmg parameters (see Mmg documentation for more information)
  hausd = Real 1000.0 !Hausdorff parameter (controls the refinement near boundaries)
  hgrad = Real 1.3  !gradation value (controls the ratio between two adjacent edges)

  verbosity = Integer 10 !Mmg verbosity

End
```
Boundary condition sections must be present to add a boundary condition Id to the new boundary elements.

```
Boundary Condition X
!! The adapted mesh boundary elements will receive the //Id//: **X**
  Target Boundaries = Y
End
```

## Example
An example for isotropic mesh adaptation can be found under [ELMER_TRUNK]/elmerice/Tests/MMG2D_Iso.
Examples for anisotropic mesh adaptation can be found under [ELMER_TRUNK]/elmerice/Tests/MMG2D_Aniso1 and [ELMER_TRUNK]/elmerice/Tests/MMG2D_Aniso2, where the mesh size is adapted using 1 or 2 variables (i.e. combining metric informations), respectively.
