# SSA inverse methods test cases

Optimisation of the basal friction coefficient.

- RUN_TWIN.sh: Run a perfect twin experiment with no noise in the input data

For the other simulations, you have to perturbe the observations 
and create the file MacAyeal_VELOCITIES_NOISE.txt in the DATA directory.

The following scripts automatically run a suite of optimisation for different values 
of the regularisation parameter and plot the L-Curve.

- RUN.sh: Use the classical regularisation that penalises first derivatives of the optimised parameter.

- RUN_TauB.sh: Use a regularisation that penalise first derivatives of the basal shear stress.

