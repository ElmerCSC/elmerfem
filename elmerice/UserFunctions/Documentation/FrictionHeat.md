# User Function getFrictionHeat
## General Information
- **USF Fortran File:** USF_GetFrictionHeating.f90
- **USF Name:** getFrictionHeat
- **Required Input Variable(s):** None

## General Description
Calculates the frictional heating at the base. Needs to be added to the geothermal heat flux. Consistent units for capacity and conductivity are needed.
Attention the keyword “friction heat” in the body section does not make reference to this heat source.
ComputeDevStressNS, ComputeNormal and a DummySolver to create the variable need to be run.

**31.1.2014: A bug was fixed in revision 6518!**
**revision 6834:** the new version with the loads should be used. User Function [getFrictionLoads](./FrictionLoads.md). The two user functions are however different since getFrictionLoads gives a result in W whereas getFrictionHeat gives a result in W/m2. Therefore, they have to be used in a different way (see the sif example).

## SIF contents
```
! needs surface normals to be provided
Solver 1
   Equation = "Normal vector"
   Variable = "Normal Vector"   
   ! in 3dimensional simulations we have 3 entries
   Variable DOFs = 3 
   !NB: does not need to actually solve a matrix
   !    hence no BW optimization needed
   Optimize Bandwidth = Logical False 
   Procedure = "ComputeNormal" "ComputeNormalSolver"
   ! if set to True, all boundary normals would be computed by default
   ComputeAll = Logical False
End
! needs Cauchy stresses
Solver 2
  Equation = String "StressSolver"
  Procedure =  File "ComputeDevStressNS" "ComputeDevStress"
!  ! this is just a dummy, hence no output is needed
!  !-----------------------------------------------------------------------
  Variable = -nooutput "Sij"
  Variable DOFs = 1
  ! the name of the variable containing the flow solution (U,V,W,Pressure)
  !-----------------------------------------------------------------------
  Flow Solver Name = String "Flow Solution"
  Exported Variable 1 = "Stress" ! [Sxx, Syy, Szz, Sxy] in 2D
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
! creates a variable that is filled with the values from the boundary condition
Solver 3
  Equation = "Dummy"
  Procedure = File "DummySolver" "DummySolver"
  Variable = String "friction"
  Variable DOFs = 1
End
Material 1
 ...
   Cauchy = Logical True !needed for ComputeDevStressNS
 ...
End
Boundary Condition 1
  Name = "bed"
  ComputeNormal = Logical True !needed for ComputeNormal
 ....
  Temp Flux BC = Logical True
  friction = Variable Coordinate 1
       Real Procedure "ElmerIceUSF" "getFrictionHeat"
  Temp Heat Flux = Variable friction
      Real MATC "0.063*(31556926.0)*1.0E-06 +tx" !assuming 63mW/m2 as geothermal heatflux
End
```

## Example
An example can be found in [ELMER_TRUNK]/elmerice/Tests/FrictionHeat.
