#Anisotropic ice rheology - AIFlow Solver

##General Information
- **Solver Fortran File:** AIFlowSolve_nlD2.f90 and AIFlowSolve_nlS2.f90
- **Solver Name:** AIFlowSolver_nlD2 and AIFlowSolve_nlS2
- **Required Output Variable(s):** AIFLow
- **Required Input Variable(s):** Temperature, Fabric
- **Optional Output Variable(s):** DeviatoricStress, StrainRate and Spin
- **Optional Input Variable(s):** None

##General Description
Solves the Stokes equation for the General Orthotropic Flow Law (GOLF) as a function of the fabric. The fabric is described using the second-order orientation tensor and its evolution can be computed using the [Fabric Solver](./FabricSolve.md). There are two different versions of the AIFlow solver depending on the non-linear extension of the flow law applied (see SIF section comments).

The anisotropic rheology as a function of the fabric is stored in a file of the type 040010010.Va. This file contains the dimensionless viscosity tabulated on a regular grid in the space spanned by the two largest eigenvectors of the second-order orientation tensor. This file is the output of a separate run of a micro-macro model (some viscosity input files can be downloaded [here](./viscosityfiles.tar.gz)). The name file (abcdefghi.Ma) contains the information about the micro-scale and type of micro-macro model used. Its nomenclature is:

- grain anisotropy parameter beta=0.abcd
- grain anisotropy parameter gamma=e.fg
- stress exponent n=h.i
- model used for tabulation =M (V holds for VPSC model)
###2.5D model – AIFlow solver accounting for flow width

Any real ensemble of flow lines may widen or get narrow, so the width of this flow tube can be accounted for in a two dimensional (x,z) model in the AIFlow solver (2.5D model). In the **Material section**, add the **FlowWidth** key word, that contains the width of the flow tube. For mass conservation, the accumulation area that should be considered correspond to the upper surface area that depends on the flow width.

##SIF contents
```
! Solve the equation for the orthotropic flow law
!  AIFlow Solvers
Solver 1
  Equation = AIFlow
  Variable = AIFlow
  Variable DOFs = 4                        !3 for 2D (u,v,p) -- 4 for 3D (u,v,w,p)


  Exported Variable 1 = Temperature        !Define Temperature Mandatory!!
  Exported Variable 1 DOFS = Integer 1

  Exported Variable 2 = Fabric             !Define Fabric Variable !!Mandatory if Isotropic=False
  Exported Variable 2 DOFS = Integer 5

  Exported Variable 3 =  StrainRate        ! Compute SR
  Exported Variable 3 DOFS = Integer 6     !4 in 2D  6 in 3D (11,22,33,12,23,31)

  Exported Variable 4 =  DeviatoricStress  ! Compute Stresses
  Exported Variable 4 DOFS = Integer 6     !4 in 2D  6 in 3D  (11,22,33,12,23,31)

  Exported Variable 4 =  Spin              ! Compute Spin
  Exported Variable 4 DOFS = Integer 3     !1 in 2D  3 in 3D (12,23,31)

! If non-linearity introduced using deviatoric stress second invariant 
  Procedure = "ElmerIceSolvers" "AIFlowSolver_nlS2"
  
! If non-linearity introduced using strain-rate second invariant   
!  Procedure = "ElmerIceSolvers" "AIFlowSolver_nlD2"
End
! Body Force
Body Force 1
  AIFlow Force 1 = Real 0.0
  AIFlow Force 1 = Real 0.0
  AIFlow Force 3 = Real -0.00899  ! body force, i.e. gravity * density
End

! Material
Material 1
!!!!! For AIFlows...
  Powerlaw Exponent = Real 3.0         ! sqrt(tr(S^2/2))^n if AIFlow_nlS2 sqrt(tr(2D^2))^(1/n-1) if  AIFlow_nlD2
  Min Second Invariant = Real 1.0e-10  ! Min value for the second invariant of strain-rates
  Reference Temperature = Real -10.0   ! T0 (Celsius)!
  Fluidity Parameter = Real 20.        ! Bn(T0)
  Limit Temperature = Real -5.0        ! TL  (Celsius)!
  Activation Energy 1 = Real 7.8e04    ! Joule/mol for T&lt;TL
  Activation Energy 2 = Real 7.8e04    ! Joule/mol for T&gt;TL

  Viscosity File = FILE "040010010.Va"

  Isotropic = Logical False   !If set to true Glen flow law (no need to define Fabric)
End

!Initial Conditions
Initial Condition 1
! Define an isotropic fabric
  Fabric 1 = Real 0.33333333333333 !a2_11
  Fabric 2 = Real 0.33333333333333 !a2_22
  Fabric 3 = Real 0.               !a2_12
  Fabric 4 = Real 0.               !a2_23
  Fabric 5 = Real 0.               !a2_13

  AIFlow 1 = Real 0.0              ! u
  AIFlow 2 = Real 0.0              ! v
  AIFlow 3 = Real 0.0              ! w
  AIFlow 4 = Real 0.0              ! p
End

! Boundary Conditions
Boundary Condition 1
  Target Boundaries = 1
!dirichlet condition for velocity
   AIFlow 1 = Real 0.0
   AIFlow 2 = Real 0.0
End


Boundary Condition 2
  Target Boundaries = 2
! Neuman condition for AIFlow
  Normal force = Real 0.0    ! force along normal
  Force 1 = Real 0.0        ! force along x
  Force 2 = Real 0.0        ! force along y
  Force 3 = Real 0.0        ! force along z

  AIFlow Slip Coeff 1 = Real 0.0   ! Slip coeff.
End
```
##Examples
[ELMER_TRUNK]/elmerice/Tests/AIFlowSolve

##References
- Extension of the linear version of the GOLF law to its non-linear form is presented in this publication:
Ma Y., O. Gagliardini, C. Ritz, F. Gillet-Chaulet, G. Durand and M. Montagnat, 2010. Enhancement factors for grounded ice and ice shelves inferred from an anisotropic ice-flow model. J. Glaciol., 56(199), p. 805-812.

- Fabric evolution and numerical implementation within Elmer/Ice are presented in this publication:
Gillet-Chaulet F., O. Gagliardini , J. Meyssonnier, T. Zwinger and J. Ruokolainen, 2006. Flow-induced anisotropy in polar ice and related ice-sheet flow modelling. J. Non-Newtonian Fluid Mech., 134, p. 33-43.

- The GOLF law is presented in detail in this publication:
Gillet-Chaulet F., O. Gagliardini , J. Meyssonnier, M. Montagnat and O. Castelnau, 2005. A user-friendly anisotropic flow law for ice-sheet modelling. J. of Glaciol., 51(172), p. 3-14.

- 2.5D model – AIFlow solver accounting for flow width:
Passalacqua O., Gagliardini O., Parrenin F., Todd J., Gillet-Chaulet F. and Ritz C. Performance and applicability of a 2.5D ice flow model in the vicinity of a dome, Geoscientific Model Development, 2016 (submitted).
