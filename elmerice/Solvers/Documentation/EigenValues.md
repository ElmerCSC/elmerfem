# Solver ComputeEigenValues
## General Information
- **Solver Fortran File:** ComputeEigenValues.f90
- **Solver Name:** ComputeEigenValues
- **Required Output Variable(s):** default is EigenStress (else in EigenValue Variable Name)
- **Required Input Variable(s):** A tensor in Tensor Variable Name
- **Optional Output Variable(s):** EigenVector1, EigenVector2 and EigenVector3
- **Optional Input Variable(s):** None

## General Description
The aim of this solver is to compute the eigenvalues of a tensor variable (e.g., strain-rate or stress tensors). Optionally, the 3 eigenvectors corresponding to each of the three eigenvalues can be computed.

## SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 3
  Equation = "EigenSR"
  Variable = -nooutput dumy
  Variable DOFs = 1

  Procedure = "ElmerIceSolvers" "ComputeEigenValues"

! 3 Eigenvalues    
  Exported Variable 1 = "EigenSR"
  EigenValue Variable Name = String "EigenSR"
  Tensor Variable Name = String "StrainRate"
  Exported Variable 1 DOFS = 3

! Principal vectors (optional) 
  Exported Variable 2 = EigenVector1
  Exported Variable 2 DOFS = 3
  Exported Variable 3 = EigenVector2
  Exported Variable 3 DOFS =  3
  Exported Variable 4 = EigenVector3
  Exported Variable 4 DOFS = 3
End
```
## Examples
An example in which the eigenvalues and eigenvectors of the strain-rate tensor are computed can be found here [ELMER_TRUNK]/elmerice/examples/Test_StrainRate.
