#Solver ComputeStrainRate
##General Information
- **Solver Fortran File:** ComputeStrainRate.f90
- **Solver Name:** ComputeStrainRate
- **Required Output Variable(s):** default is StrainRate (else in StrainRate Variable Name)
- **Required Input Variable(s):** A Flow Solution (in Flow Solution Name)
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

##General Description
The aim of this solver is to compute the strain-rate tensor from the flow solution. For a 2D simulation there are 5 DOFs (E11, E22, E33, E12, Eii) and for a 3D simulation, 2 additional (7 DOFs) are being solved for (E11, E22, E33, E12, E23, E31, Eii). This solver uses a dummy variable and solves 5 (7 in 3D) times a 1 DOF system for each strain-rate components and the trace of the strain-rate tensor (Eii).

The Strain-rate tensor is computed using:
*epsilon_{ij}  =  {1/2} ({u_{i,j} + u_{j,i}})*
where *u* is the velocity vector solution of the Stokes problem.

##SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 1
  Equation = "Strain Rate"
  Procedure = "ElmerIceSolvers" "ComputeStrainRate"
! this is just a dummy, hence no output is needed
!-----------------------------------------------------------------------  
  Variable = -nooutput "Eij"
  Variable DOFs = 1

  Exported Variable 1 = "StrainRate"
  Exported Variable 1 DOFs = 5 !in 2D, 7 in 3D

! the name of the variable containing the flow solution (U,V,W,Pressure)
!-----------------------------------------------------------------------
  Flow Solver Name = String "Flow Solution"
! the name of the strain-rate solution (default is 'StrainRate')
  StrainRate Variable Name = String "StrainRate"
  
  Linear System Solver = Direct
  Linear System Direct Method = umfpack
End
```

##Examples
An example using this solver can be found in [ELMER_TRUNK]/elmerice/examples/Test_StrainRate.
