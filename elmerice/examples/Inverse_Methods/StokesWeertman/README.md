# Inverse methods test cases
- Creation date: 28/02/2023
- Modification date : 28/02/2023

## Description  

Example for the optimisation of the Weertman friction coefficient using the vectorised incompressible Stokes solver (**IncompressibleNSVec**).

This is a twin experiment where we try to recover the true friction coeffcient, using the solution computed with the model as target.


## Run the example

1. Make mesh 

	```
	ElmerGrid 1 2 rectangle
	```

2. Generate solution 

	```
        ElmerSolver Direct_nl.sif
        ```

3. Run the optimisation

	```
        ElmerSolver OPTIM_TWIND_nl.sif 
        ```


