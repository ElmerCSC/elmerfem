# Solver ComputeDevStress
## General Information
- **Solver Fortran File:** ComputeDevStressNS.f90
- **Solver Name:** ComputeDevStress
- **Required Output Variable(s):** default is Stress (else in Stress Variable Name)
- **Required Input Variable(s):** A Flow Solution (in Flow Solution Name)
- **Optional Output Variable(s):** None
- **Optional Input Variable(s):** None

## General Description
The aim of this solver is to compute deviatoric or Cauchy stress from flow solution. For a 2D simulation there are 4 DOFs (S11, S22, S33, S12), for a 3D simulation, 2 additional are being solved for (S11, S22, S33, S12, S23, S31). This solver uses a dummy variable and solves 4 (6 in 3D) times a 1 DOF system for each stress components.

The Cauchy stress is computed using:
*sigma_{ij}  = 2 {eta}  {epsilon}_{ij} - p delta_{ij}*
where *epsilon* is directly evaluated from the velocity field and *p* is the isotropic pressure.The convention is that a positive stress corresponds to a tensile stress (opposite to the isotropic pressure convention).

This solver doesn't work for the GOLF anisotropic ([AIFlow Solver](./AIFlowSolve.md)) rheology. Nevertheless, this solver has intrinsic functions that allow direct computation of the stress.

2024-01-23: the ComputeDevStress solver is now compatible with the snow/firn rheology ([Porous Solver](./PorousSolve.md)). In that case, the Flow Solver Name variable has simply to be set to "Porous". Note that the possibility of using intrinsic functions of the Porous Solver allowing direct computation of stresses at nodes is kept. In that case, the stress at a given node is computed as the average contribution from all the elements belonging at this given node. This is slightly different than solving the above equation through the variational method as done by the ComputeDevStress. The ComputeDevStress approach is preferable, especially at partition interfaces.

2023-11-14: Using the new Vectorized Stokes solver IncompressibleNSvec, the ComputeDevStress solver will not work as the Constant Temperature in the Material has been changed to Relative Temperature. You can either set back to the old keywords setting Glen Allow Old Keywords to True in the Material section. An other solution is to use the effective viscosity and strain-rate computed internaly in the vectorized Stokes solver (variables Viscosity and Shearrate should be exported using either -ip or -lm).  

## SIF contents
The required keywords in the SIF file for this solver are:

```
Solver 1
  Equation = String "StressSolver"
  Procedure =  File "ElmerIceSolvers" "ComputeDevStress"
  ! this is just a dummy, hence no output is needed
  !-----------------------------------------------------------------------
  Variable = -nooutput "Sij"
  Variable DOFs = 1
  ! the name of the variable containing the flow solution (U,V,W,Pressure)
  !-----------------------------------------------------------------------
  Flow Solver Name = String "Flow Solution" ! "Porous" is if used together with Porous Solver
  ! no default value anymore for "Stress Variable Name"
  Stress Variable Name = String 'Sigma'
  !-----------------------------------------------------------------------
  Exported Variable 1 = "Sigma" ! [Sxx, Syy, Szz, Sxy] in 2D
                                 ! [Sxx, Syy, Szz, Sxy, Syz, Szx] in 3D
  Exported Variable 1 DOFs = 6   ! 4 in 2D, 6 in 3D
  Linear System Solver = "Iterative"
  Linear System Iterative Method = "BiCGStab"
  Linear System Max Iterations = 300
  Linear System Convergence Tolerance = 1.0E-09
  Linear System Abort Not Converged = True
  Linear System Preconditioning = "ILU0"
  Linear System Residual Output = 1
End

Material 1
  ...
  ! we want to have the Cauchy stress
  !----------------------------------
  Cauchy = Logical True
End
```

## Examples
A 2D example can be found in [ELMER_TRUNK]/elmerice/Tests/ComputeDevStress.

## Reference
This solver can be cited using the following references:
Gagliardini O., D. Cohen, P. Råback and T. Zwinger, 2007. Finite-Element Modeling of Subglacial Cavities and Related Friction Law. J. of Geophys. Res., Earth Surface, 112, F02027.
