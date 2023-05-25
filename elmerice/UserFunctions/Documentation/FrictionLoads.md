# User Function getFrictionLoads
## General Information
- **USF Fortran File:** USF_GetFrictionHeating.f90
- **USF Name:** getFrictionLoads
- **Required Input Variable(s):** None

## General Description
The function getFrictionLoads replaces the old function getFrictionHeat, which still remains in the svn because of compatibility issues (from **revision 6834** on). It is recommended to switch to this newer version. The two user functions are however different since getFrictionLoads gives a result in W whereas getFrictionHeat gives a result in W/m2. Therefore, they have to be used in a different way (see the sif example).
This approach is based on the coupling of surface heat production to the residual of the Stokes solver - the first three components of the residual represent the nodal forces acting on the surface. This is the natural way of coupling in FEM.
Note, in the case of a restart-run, the (Navier-)Stokes Solver needs to exist in the sif-file (option Exec Solver = Never is sufficient), otherwise the Loads are not found.
Note, the function has not yet been tested in the case of basal melting, but should be working.
Sometimes negative values for the friction loads are calculated and set to zero in the source code. This has not yet been fully understood and could lead to some instabilities. This has, however, not yet been observed.

Output can be produced by using the ForceToStress Solver (see below).

Consistent units for capacity and conductivity are needed.
Attention, the keyword “friction heat” in the body section does not make reference to this heat source - but to “strain heat” or “deformational heat”.

## SIF contents
```
!Normals are needed 
Solver x 
  Equation = "NormalVector" 
  Exec Solver = Before Simulation
  Procedure = "ElmerIceSolvers" "ComputeNormalSolver"
  Variable = String "Normal Vector" 
  Variable DOFs = 2
  Optimize Bandwidth = Logical False 
  ComputeAll = Logical False
End

!add in NS Solver Section
  Exported Variable 1 = Flow Solution Loads[Fx:1 Fy:1 Fz:1 CEQ Residual:1 ]
  Calculate Loads = Logical True

!add to Boundary Conditions BED
   Mass Consistent Normals = Logical True 
   ComputeNormal = Logical True !needed for ComputeNormal
   
   Temp Load = Variable Velocity 1
          Real Procedure  "ElmerIceUSF" "getFrictionLoads"	

   !geothermal heat flux (as usual - not modified as it was the case with the earlier version)
   Temp Flux BC = Logical True (needed when geothermal flux is used)
   Temp Heat Flux = Variable Coordinate 1
        Real MATC "0.040*(31556926.0)*1.0E-06" !40mW/m2
```
## Output of the friction heat
Output can be produced with the ForceToStress Solver. Sometimes negative values are calculated, which should be set to zero since they are only numerical artefacts. Bear in mind that the variable names here are just for output and can be changed, but should not interfere with naming of the Temperature load variable.

```
!variable created for output purpose, to be executed at the beginning ONLY on the BED
Solver x 
  Equation = "Dummy"
  Procedure = File "DummySolver" "DummySolver"
  Variable = String "Friction Load" !variable created for output purpose
  Variable DOFs = 1
End

!after NS, to be executed ONLY on the BED
Solver x 
  Equation = "ForceToStress"
  Procedure = File "ElmerIceSolvers" "ForceToStress"
  Variable = String "Friction Heating"
  Variable DOFs = 1
  Force Variable Name = String "Friction Load"
  Linear System Solver = "Iterative"
  Linear System Iterative Method = "BiCGStab"
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0E-05
  Linear System Abort Not Converged = False
  Linear System Preconditioning = "ILU0"
  Linear System Residual Output = 1
End

!add to Boundary Conditions BED
 Friction Load = Variable Velocity 1            
         Real Procedure  "ElmerIceUSF" "getFrictionLoads"
```

## Example
An example can be found in [ELMER_TRUNK]/elmerice/Tests/FrictionHeat.
