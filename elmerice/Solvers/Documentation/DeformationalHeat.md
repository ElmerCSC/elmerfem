#DeformationalHeat Solver
##General Information
- **Solver Fortran File:** DeformationalHeat.f90
- **Solver Name:** DeformationalHeatSolver
- **Required Output Variable(s):** W (User Defined)
- **Required Input Variable(s):** Default Flow Solution (else in Flow Solver Name)
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

##General Description
This solver computes the volumetric heat produced by ice flow deformation.

##SIF contents
This solver produces the volumetric heat source due to ice deformation (strain heating). This can be directly used within the TemperateIceSolver using the body-force keyword Temp Volume Source (provided your temperature variable is called Temp. Important: Do not use that as a Heat Source for the regular HeatSolve (it would need a specific heat source).

The required keywords in the SIF file for this solver are:

```
Solver 3
  Equation = DeformationalHeat
  Variable = W
  Variable DOFs = 1

  procedure =  "ElmerIceSolvers" "DeformationalHeatSolver"

  Linear System Solver = direct
  Linear System direct Method = umfpack
  
  Flow Solver Name = String 'Flow Solution'
End

! This would be the correct source for TemperateIceSolver
! with variable name Temp
Body Force 1
   Temp Volume Source = Equals W
End
```
##Example
A test using the DeformationalHeat solver can be found in [ELMER_TRUNK]/elmerice/Tests/Teterousse_DeformHeat.
