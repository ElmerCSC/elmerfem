# User Function USF_ShapeFactor
## General Information
- **USF Fortran File:** USF_ShapeFactor.f90
- **USF Name:** LateralFriction_x and LateralFriction_y
- **Required Input Variable(s):** None

## General Description
This user function return the modified gravity force *-g +(1-f)(g.t).t*
where *t* is the tangent vector to the upper surface and *f* the shape factor. The shape factor *f* can be given or calculated as a function of *h(x,t)* for rectangular or parabola valley shapes.

This solver works only in 2D (no sense in 3D). It works for non-structured mesh.

## SIF contents
Here are the three options to use the USF_ShapeFactor user function. In the first option, the shape factor *f* is simply given as a parameter Shape Factor:

```
Body Force 1
  Flow BodyForce 1 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "ShapeFactorGravity_x"
  Flow BodyForce 2 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "ShapeFactorGravity_y"
  Shape Gravity 1 = Real  0.0                                 
  Shape Gravity 2 = Real -9.7696e15       

  Shape Factor = real $shapefactor
End
```
In the second option, the shape factor *f* is computed for a rectangular valley of given width:

```
Body Force 1
  Flow BodyForce 1 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "ShapeFactorGravity_x"
  Flow BodyForce 2 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "ShapeFactorGravity_y"
  Shape Gravity 1 = Real  0.0                                 
  Shape Gravity 2 = Real -9.7696e15       

  Rectangle Shape = Logical True
  Rectangle Width = Real 2.0e3
End
```
In the third option, the shape factor *f* is computed for the parabola type transverse valley *y = b + a*z^2*, where *a* is given by the Parabola afactor keyword and *b* is the bedrock elevation of the central line:

```
Body Force 1
  Flow BodyForce 1 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "ShapeFactorGravity_x"
  Flow BodyForce 2 = Variable Coordinate 1
     Real Procedure "ElmerIceUSF" "ShapeFactorGravity_y"
  Shape Gravity 1 = Real  0.0                                 
  Shape Gravity 2 = Real -9.7696e15       
  
  Parabola Shape = Logical True
  Parabola afactor = Real 0.0005
End
```
In all cases, the two keywords Shape Bedrock and Shape Surface are used to locate nodes belonging on the bedrock and surface boundaries, respectively.

```
! Bedrock
Boundary Condition 1
  Target Boundaries = 1
  ...
  Shape Bedrock = Logical True
End

! Upper Surface
Boundary Condition 3
  Target Boundaries = 3
  ...
  Shape Surface = Logical True
End
```

## Examples
An example using the user function USF_ShapeFactor can be found in [ELMER_TRUNK]/elmerice/Tests/ShapeFactor.

## Reference
When used this solver can be cited using the following references:
Jay-Allemand M., F. Gillet-Chaulet, O. Gagliardini and M. Nodet, 2011. Investigating changes in basal conditions of Variegated Glacier prior to and during its 1982–1983 surge. The Cryosphere, 5, p. 659-672, doi:10.5194/tc-5-659-2011.
