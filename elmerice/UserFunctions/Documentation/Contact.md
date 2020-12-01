# User Function USF_Contact
## General Information
- **USF Fortran File:** USF_Contact.f90
- **USF Name:** SlidCoef_Contact
- **Required Input Variable(s):** GroundedMask, Fw, Flow Solution Loads, Normal Vector
- **Optional Input Variable(s):** Distance

## General Description
The aim of this user function is to test the contact during the non-linear iteration of the Stokes solver. This is done through the application of basal sliding. If the Grounding line is retreating, nodes can move from grounded to floating and the variable GroundedMask is updated accordingly. The contact is tested once after the solution of the non-linear iterations has slightly converged and then is kept fixed until the complete convergence of the non-linear iterations.

This user function uses the variables GroundedMask, computed by [GroundedSolver](../../Solvers/Documentation/GroundedSolver.md), Fw, computed by [GetHydrostaticLoads](../../Solvers/Documentation/GetHydrostaticLoads.md), the residual of the Stokes solver in Flow Solution Loads, and the Normal Vector from the [ComputeNormal](../../Solvers/Documentation/Normal.md) Solver. If the keyword Non Detachment Inland Distance has a positive value, the variable Distance computed by the Distance Solver of the Elmer distribution is also needed.

## SIF contents
The required keywords in the SIF file for these user functions are:

```
!! BC  Bedrock + Shelf
Boundary Condition 1
  Name = "bottom"
  Target Boundaries = 1
  Body Id = 3
  
  ...
  Grounding Line Moves = Logical True ! Default is true, 
                                      !useful to test influence of not moving the GL
  
  Non Detachment Inland Distance = 10000.0 ! distance from the GL where nodes 
                                           ! are fixed to the bed for the Free Surface evolution
                                           ! Default is -1000.0 (no nodes fixed) 
                                           ! If > 0.0, needs the Distance Solver
  
  Slip Coefficient 3 = Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "SlidCoef_Contact"
  Slip Coefficient 2 = Variable Coordinate 1
    Real Procedure "ElmerIceUSF" "SlidCoef_Contact"

  Sliding Law = String "Weertman"  ! Alternative is "Coulomb"
  Weertman Friction Coefficient = Real $C
  Weertman Exponent = Real $(1.0/n)
  Weertman Linear Velocity = Real 1.0
  ! Options are 'Last Grounded' (default), 'First Floating' or 'Discontinuous' 
  Grounding Line Definition = String "Discontinuous"
  Test Contact Tolerance = real 1.0e-3
End
```
Depending on the chosen sliding law, the corresponding parameters must be given. See the documentation for the corresponding user functions Weertman, Budd and Coulomb.

## Examples
2D examples can be found in [ELMER_TRUNK]/elmerice/Tests/Contact and [ELMER_TRUNK]/elmerice/Tests/GL_MISMIP.
