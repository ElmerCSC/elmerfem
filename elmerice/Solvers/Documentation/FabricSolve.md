# Fabric Evolution - DG Fabric Solver
## General Information
- **Solver Fortran File:** FabricSolve.f90
- **Solver Name:** FabricSolver
- **Required Output Variable(s):** Fabric
- **Required Input Variable(s):** Temperature, AIFlow,
- **Optional Output Variable(s):** EigenV
- **Optional Input Variable(s):** Mesh Velocity

## General Description
This solver solves the fabric evolution equations for the anisotropic law obtained from the [AIFlow solver](./AIFlowSolve.md). In Elmer/Ice, the fabric is described using the 5 independent components of the second-order orientation tensor. The grain rotation is induced by the macroscopic strain-rate and the deviatoric stress, and their relative influence is controlled by the interaction parameter (from 0 for purely strain-induced rotation to 1 for purely stress induced rotation). To account for other phenomena than grain rotation, one can use the diffusion parameter.

The eigenvalues of the second-order orientation tensor Fabric are computed if the variable EigenV exists.

## SIF contents
```
! Solve the equation for the orthotropic flow law
!  AIFlow Solvers
Solver 1
  Equation = AIFlow
  Variable = AIFlow
  Variable DOFs = 4                        !3 for 2D (u,v,p) -- 4 for 3D (u,v,w,p)


  Exported Variable 1 = Temperature        ! Define Temperature Mandatory!
  Exported Variable 1 DOFS = Integer 1

  Exported Variable 2 = Fabric             ! Define Fabric Variable
  Exported Variable 2 DOFS = Integer 5     ! Mandatory if Isotropic=False


! If non-linearity introduced using deviatoric stress second invariant 
  Procedure = "ElmerIceSolvers" "AIFlowSolver_nlS2"
End

! Fabric solver itself
Solver 2
  Equation = Fabric
  Variable = -nooutput Compfab    ! dumy variable
  Variable DOFs = 1               !FabricSolver compute each variable independently, 
                                  !Picard Type iterations

  Procedure = "ElmerIceSolvers" "FabricSolver"
  Discontinuous Galerkin = Logical True
End
! Material
Material 1
!!!! For Fabric Solver
  Interaction Parameter = Real 0. ! 0 => Fabric Evolution function of Strain-rates 
                                  ! 1 => Fabric Evolution function of dev stresses
                                  !If not defined set to the default value given in the Viscosity File
                                  
  Diffusion Parameter = Real 0.   ! Diffusion term. By default set to 0 if not defined
End

!Initial Conditions
Initial Condition 1
! Define an isotropic fabric
  Fabric 1 = Real 0.33333333333333 !a2_11
  Fabric 2 = Real 0.33333333333333 !a2_22
  Fabric 3 = Real 0.               !a2_12
  Fabric 4 = Real 0.               !a2_23
  Fabric 5 = Real 0.               !a2_13
  ...
End

! Boundary Conditions
Boundary Condition 1
  Target Boundaries = 1
! Dirichlet conditions for Fabric 
! only if inflow boundary condition, no condition for outflow
  Fabric 1 = Real 0.33333333333333 !a2_11
  Fabric 2 = Real 0.33333333333333 !a2_22
  Fabric 3 = Real 0.               !a2_12
  Fabric 4 = Real 0.               !a2_23
  Fabric 5 = Real 0.               !a2_13
End
```

## Examples
Download an example using the Fabric Solver. TODO

## References
- Fabric evolution and numerical implementation within Elmer/Ice are presented in this publication:
Gillet-Chaulet F., O. Gagliardini , J. Meyssonnier, T. Zwinger and J. Ruokolainen, 2006. Flow-induced anisotropy in polar ice and related ice-sheet flow modelling. J. Non-Newtonian Fluid Mech., 134, p. 33-43.
