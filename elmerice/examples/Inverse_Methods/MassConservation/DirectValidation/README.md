# Direct Validation 

Validate the model results against the analytical solutions.
	- The SSA is solved using the analytical ice thickness 
	- The steady state thickness is computed using the analytical velocity and SMB

SaveScalars save the rms of h-h^{true} and u-u^{true}; the error should decrease as the resolution increases.  

We test the convergence using linear and quadratic elements

## Content

* *MakeValidation.sh*: bash script that generate meshes for different resolution and run the simulations 
* *RAMP.sif*  *RAMP_2nd.sif*: .sif files to compute the solution


