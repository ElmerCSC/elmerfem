## Icy Mask Solver  {#IcyMaskSolver}

**Module name**: IcyMaskSolver.F90<br>
**Module subroutines**: IcyMaskSolver<br>
**Module authors**: Olivier Gagliardini (IGE-Grenoble)<br>
**Document authors**: Olivier Gagliardini<br> 
**Document edited**: 06/12/2020<br>

**Required input variables:**

 - Thickness (Computed using #StructuredProjectToPlane   

**Output variables:**

 - IcyMask 

**Keywords:**

 - 'Toler' [Real] tolerance for testing (Thickness - Toler > Hmin) [set in Solver]
 - 'Ice Free Thickness' (Real) - Hmin [set in Solver]
 - 'Remove Isolated Points' (Logical) - to impose a Mask < 0 to isolated nodes surounded by nodes with all Mask <0 [set in Solver]
 - 'Remove Isolated Edges' (Logical) - to impose a Mask < 0 to isolated edges surounded by nodes with all Mask <0 [set in Solver]

### Introduction

This solver is to be used to locate Ice Free area of a domain. Ice Free areas are defined by area where Thickness is lower than a minimum thickness. 

The Mask value is 
 - +1 for glacier nodes (H > Hmin)
 - 0 for the nodes belonging to the glacier contour (H=Hmin)
 - -1 for Ice free nodes (H <= Hmin)
 - <-1 for isolated nodes 

Isolated nodes (and edges) are defined by nodes (edges) surounded by only ice free nodes. They are tagued with a mask value < -1.  

#### SIF 

```
Solver 3
  Equation = "HeightDepth"
  Procedure = "StructuredProjectToPlane" "StructuredProjectToPlane"
  Active Coordinate = Integer 3
  Operator 1 = Thickness
End

Solver 4 
  ! to be executed on top surface (need Thickness)
  Equation = "IcyMask"
  Variable = "IcyMask"
  Variable DOFs = 1
  Procedure = "bin/IcyMaskSolver" "IcyMaskSolver"

  Toler = Real 1.0e-1
  Ice Free Thickness = Real #MinH
  Remove Isolated Points = Logical True 
  Remove Isolated Edges = Logical True ! Need that Remove Isolated Points is set to True to be accounted for 
End 
```

