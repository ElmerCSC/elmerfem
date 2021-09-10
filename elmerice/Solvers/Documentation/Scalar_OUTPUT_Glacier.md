## Scalar\_OUTPUT\_Glacier Solver  {#Scalar\_OUTPUT\_Glacier}

**Module name**: Scalar\_OUTPUT\_Glacier.F90  
**Module subroutines**: Scalar\_OUTPUT  
**Module authors**: Olivier Gagliardini (IGE-Grenoble)  
**Document authors**: Olivier Gagliardini  
**Document edited**: 
	- 21/12/2020:  
	- 04/05/2021: F. Gillet-Chaulet  
		- Integrations now done on elements projected on the horizontal plane
		- Add output of volume change (integrate surface elevation rate of change) and residual ice flux due to limits and Dirichlet conditions.

**Required input variables:**

 - Thickness (Computed using #StructuredProjectToPlane)   
 - Zs (From the free surface solver)
 - BedDEM (basal surface elevation)
 - Zs Velocity (surface elevation rate of change; computed by adding *Calculate Velocity = Logical True* in the FreeSurface solver)
 - Zs Loads (residual due to limits and Dirichlet conditions; N.B. Zs residual only outputs residual due to enforced limits; computed using *calculate loads = Logical True*  in the FreeSurface solver)

**Required Material Properties:**  
  - *Min Zs*: the lower limit imposed to the FreeSurface solver.  

**Required Body Forces:**  
  - *Zs Accumulation Flux 3*: the surface mass balance forcing of the FreeSurface solver.



**Output variables:**

 - None, output scalar variable in an ascii file

**Keywords:**

 - 'File Name' [String] the name of the ascii output file [set in Solver]
 - 'XRefFrontPosition' and 'YRefFrontPosition' [Real] the initial front position X and Y coordinates [set in Solver]
  

### Introduction

This solver is to be used to output some scalar quantities for a glacier configuration (domain without ice characterised by an IcyMask < 0). 
The quantities are: glacier volume, glacier area, Ablation area, Accumulation area, SMB total, SMB Ablation, SMB Accumulation and Front elevation. 

#### SIF 

```
Solver ....
  Equation =  String "Free Surface Evolution"
  Procedure = "FreeSurfaceSolver" "FreeSurfaceSolver"

  Variable = "Zs"
  Variable DOFs = 1

  ! Apply internal limiters
  Apply Dirichlet = Logical true

  ! calculate dz/dt (better than from mesh velocity in case of steady-state iterations)
  Calculate Velocity = Logical True

  ! loads also takes into account dirichlet conditions
 ! to compute residual flux
  calculate loads = Logical True
End 
```


```
Solver 10
  Exec Solver = After Timestep
  Equation = "Save 1D Vars"
  Procedure = File "bin/Scalar_OUTPUT_Glacier" "Scalar_OUTPUT"
  Variable = -nooutput "savescal"

  File Name = File "1DVar_OUTPUT_$namerun$.dat"

  XRefFrontPosition = real 9.5635e+05
  YRefFrontPosition = real 1.1916e+05
End
```

