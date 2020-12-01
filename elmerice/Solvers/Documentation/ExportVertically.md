#Solver ExportVertically
##General Information
- **Solver Fortran File:** ExportVertically.f90
- **Solver Name:** ExportVertically
- **Required Output Variable(s):** user defined
- **Required Input Variable(s):** None
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

##General Description
The aim of this solver is to export vertically in the whole domain the value of the variable set as a Dirichlet condition on the upper or lower boundary. This can be used for example if one want to know in the whole domain the GroundedMask variable which is only computed on the bedrock boundary.

##SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 1
  Equation = "ExportVertically"
  Procedure = File "ElmerIceSolvers" "ExportVertically"
  Variable = String "ExportedMask"
  Variable DOFs = 1
  Linear System Solver = "Direct"
  Linear System Direct Method = umfpack
End

Boundary Condition 1
  Target Boundaries = 1
  ExportedMask = Equals GroundedMask 
End
```

##Examples
An example in which an analytical function is exported vertically on the whole domain can be found in [ELMER_TRUNK]/elmerice/Tests/ExportVertically.
