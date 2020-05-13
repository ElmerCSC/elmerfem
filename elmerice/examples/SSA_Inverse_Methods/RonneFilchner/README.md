# Ronne-Filchner these case

Optimisation of the viscosity for the Ronne-Filchner ice shelf.
To run this experiment you first need to download observations of the surface velocities, see instructions in the DATA directory.

- RUN.sh: run a suite of optimisations with different values
of the regularisation parameter and plot the L-Curve.

- mesh: contain the mesh and elmer .result files (IMPORT.result.*) where the variables required to run the simulation
(topography, initial viscosity filed) have been initialised.

