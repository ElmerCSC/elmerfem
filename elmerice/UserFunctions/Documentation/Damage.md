# User Function USF_Damage
## General Information
- **User Function Fortran File:** USF_Damage.f90
- **Function 1 Name:** EnhancementFactor
- **Function 2 Name:** SourceDamage

## General Description
[USF_Damage.pdf](./usf_damage.pdf)

The first function changes the rheology of the damage ice by impacting the enhancement factor.

The second function is used to compute the source term of the advection-reaction solver.

## SIF contents
The required keywords in the SIF file for this User Function are:

```
Constants
  Dev Tensile Strength Modifier = Real 0.05 ! standard deviation for the stress threshold distribution
End

Body Force 1
  DGD Source = Variable Damage
    Real Procedure "ElmerIceUSF" "SourceDamage"
End

Material 1
  Glen Enhancement Factor = Variable Damage
    Real Procedure "ElmerIceUSF" "EnhancementFactor" 
    
  Damage Enhancement Factor = Real 2.00 ! damage enhancement factor
  Damage Parameter sigmath = Real 0.05 ! stress threshold for damage increase
End
Additionally, for output visualisation, the damage criterion Chi is saved as a variable named Chi, which need to be exported in a solver, such as :

Solver 3
  Equation = Sij
  Procedure = "ElmerIceSolvers" "ComputeDevStress"
[...]
  Exported Variable 1 = Stress[Sxx:1 Syy:1 Szz:1 Sxy:1 Syz:1 Sxz:1]
  Exported Variable 1 DOFs = 6
  Exported Variable 2 = -dofs 1 "Chi"
[...]
End
```

## Examples
A 3D example can be found in [ELMER_TRUNK]/elmerice/Tests/Damage.

## Reference
Krug, J., J. Weiss, O. Gagliardini and G. Durand, 2014. Combining damage and fracture mechanics to model calving, The Cryosphere, 8, 2101-2117, doi:10.5194/tc-8-2101-2014.
