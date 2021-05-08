# ISCAL (Ice Sheet Coupled Approximation Levels)
ISCAL couples the SIA with the full Stokes (FS) equations, so that the FS equations are used only where the SIA is inaccurate. The SIA error is estimated automatically. The core of the code is FlowSolveSIAFS.f90 which couples the SIA and FS solution, the SIASolverOnExtruded.f90 which computes the SIA-solution, and the ErrorEstimationSubs.f90 which estimates the error.

## Coupling
The coupling simply consists of setting the SIA solution as a boundary condition to the FS solution. The SIA solution is computed over the whole domain (computational cost neglible). Technically, this is done by splitting the FEM-matrix A into a block-matrix [A_FS, A_Co; A_Oc, A_SIA], where A_FS represents the Stokes equations in the FS areas, and A_SIA represents the Stokes equations in the SIA areas (i.e., A_SIA is not needed), A_Co representes the coupling from SIA to FS areas and A_Oc represents the coupling from FS to SIA areas. Only A_FS and A_Co is assembled. A_Co is multiplied by the forcing and moved to the right hand side. The remaining system A_FS * FSsolution = RHS is solved. Finally, the complete solution is constructed by combining the FSsolution in the FS areas with the SIA solution in the SIA areas. All of this is done from within FlowSolveSIAFS.f90, which is a modified version of the Elmer core routine FlowSolve.f90.

## Error Estimators
The FS is used where the SIA error is too high. There are three ways of estimating the SIA error:
- Computing the relative error in the SIA velocity.
- Computing the error in a functional of the SIA solution (such as flux through a certain point). Only implemented in 2D.
- Computing the error in the residual, i.e. how well the SIA solution fulfills the FS equations.

This is done in ErrorEstimationSubs.f90

## SIA solution
The SIASolverOnExtruded.f90 solves the SIA equations on extruded meshes in an efficient way. It is possible to instead use the built-in SIA solver or any other approximative solver (as long as it is computationally inexpensive).

## Use
The code is compiled with Compile.sh (change the path where it puts the binaries…). The code for the ISCAL is used in the same way as the rest of Elmer, i.e. you add solvers in your sif-file. The FlowSolveSIAFS.f90 is used in the same way as the original FlowSolve.f90 routine, but with some extra keywords. An example sif-entry is found below. The SIA-solver requires a few solvers to run before it (namely several FluxSolvers to compute gradients of ice surface etc.). An example of sif-entries for this is found below. The SIA-solver should always run before the coupler (FlowSolveSIAFS.f90).

There are two extra output variables for the coupling code which you can view in, for example, Paraview - the ApproximationLevel which is 1 for SIA areas and 2 for full Stokes areas, and the SIAerror, which is the error in the SIA (can be computed in three different ways).

### SIA Solver keywords description
*Velocity Cutoff = Real 50000 Sets a maximum limit on the velocity, in this case 50 000 meters per year

*Active Coordinate = 3 The coordinate direction in which the integrals of SIA will be computed - always set to the vertical direction!

*Surface Name = String FS The name of the free surface

*Bedrock Data = Logical True True if you have some data for the bedrock, false if the bedrock is flat

*Bedrock Name = String BedrockElevation The name of your variable containing the bedrock topography

!in case you use the thicknesssolver instead of free surface solver: Compute Surface From Thickness = Logical True

*Thickness Name = String H

### Coupler Keywords
- Couple Approximations = Logical True True if you wanna couple SIA and full Stokes, false if you wanna run only full Stokes

- SIA as initial condition = Logical True True if you want to use SIA as initial condition in the first timestep

- Error Estimation Intervals = Integer 20 How often you want to do error estimation, in this case it will be every 20th timestep

- Error Estimation Method = String “solution” The options are in theory “solution”, “functional” and “residual” It will estimate the error based on the horizontal velocity error directly, in a fucntional of the solution or in the residual of the solution.

**Keywords for solution based estimate**
- Relative Error Allowed In Percent = Real 15.0 Maximum allowed relative error of SIA in percent, in this case 15 %

- Absolute Error Allowed = Real 1.0 !m/a Maximum allowed absolute error in SIA in meters per year, in this case 1 meter per year

**Keywords for functional based error estimate (code will be updated)**
- Functional = String “flux across point” !The functional you are interested in. Currently only this one is implemented but you can make one yourself

- Point x-coord = Real 500000 ! For “flux across point” this is the point over which the flux is computed

- Active Coordinate = 2 ! The direction in which the flux integration is done.

- Nodewise limit for dual problem = Real 10.0 ! The allowed value of x^T*residual in each node point, where x is the dual solution

**Keywords for residual based error estimation (code will be updated)**
- Maximum Allowed Residual = Real 10000.0

**In Materials section**
If you want an isothermal ice sheet, you have to use the power law viscosity model.

- Isothermal = Logical True !read in SIA solver

For a non isothermal case, you need Arrhenius Factors

- Isothermal = Logical False

- dArrheniusFactordT = Real 0.0 ! Read in SIA Solver

This entry is used in the Navier Stokes solver on timestep 1, to say whether or not a previous estimate of the SIA error has been calculated and should be used.

- Previous ApproximationLevel = Logical False

in the case of false, the program reads this in as a dummy value for the SIA approximation level (if values are greater than zero, then it assumes the SIA solution at the node is ok)

- SIA node = Real -1.0 !

**Boundary conditions**
It is possible to put sliding conditions also for the SIA.
