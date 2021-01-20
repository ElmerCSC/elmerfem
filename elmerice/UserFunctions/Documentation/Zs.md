# User Function USF_Zs
## General Information
- **USF Fortran File:** USF_Zs.f90
- **USF Name:** ZsIni, ZsMZsIni, ZsTopIni, ZsTopMZsIni,…
- **Required Input Variable(s):** Zs or Zs Top or Zs Bottom

## General Description
The aim of these user functions is to get the appropriate value of the mesh update from the surface displacement.

The first type of user function (~Ini) is used to initialize the free surface elevation (surface elevation equals to the altitude of the nodes belonging on this surface). The second type (~MZsIni) calculate the mesh update variable as the new surface elevation minus the initial elevation.

As presented here, these user functions can be easily replaced by MATC functions (recalling that a MATC function consumes more cpu than a user function).

## SIF contents
The required keywords in the SIF file for this solver are:

```
Initial Condition 2
  Zs Top = Variable Coordinate 3
    Real Procedure "ElmerIceUSF" "ZsTopIni"
End

Initial Condition 3
  Zs Bottom = Variable Coordinate 3
    Real Procedure "ElmerIceUSF" "ZsBottomIni"
End

! Bottom boundary condition
Boundary Condition 1
  Target Boundaries = 1
!!! this BC is equal to body no. 3 !!!
  Body Id = 3
  
  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
  Mesh Update 3 = Variable Zs Bottom
    Real Procedure "ElmerIceUSF" "ZsBottomMzsIni"
End

! Upper Surface
Boundary Condition 2
  Target Boundaries = 2
!!! this BC is equal to body no. 2 !!!
  Body Id = 2

  Mesh Update 1 = Real 0.0
  Mesh Update 2 = Real 0.0
  Mesh Update 3 = Variable Zs Top
    Real Procedure "ElmerIceUSF" "ZsTopMzsIni"
End
```

## Examples
Examples of the usage of these user functions can be found in the [Elmer/Ice course material](http://elmerfem.org/elmerice/wiki/doku.php?id=courses:courses) (ISMIP application for one upper free surface and Tête Rousse application for an upper and lower free surfaces).
