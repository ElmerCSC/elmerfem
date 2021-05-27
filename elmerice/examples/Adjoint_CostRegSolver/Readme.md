# Validation test case for Adjoint_CostRegSolver:

## Available test cases
- Valid_CostRegSolver1.sif
 Test penalising first spatial derivatives of the optimised variable

- Valid_CostRegSolver2.sif
 Same but testing a change of variable and nodal value provided from the body force

- Valid_CostRegSolver3.sif
 Same but solver executed on the bottom surface of a 3D mesh

- Valid_CostRegSolver4.sif
 Regularistaion from **prior** value

- Valid_CostRegSolver5.sif
 Optimisation using a **prior**. 

## Running the tests

All cases can be run in serial or parallel

```bash
# Make mesh
## serial
ElmerGrid 1 2 mesh2D
## parallel
ElmerGrid 1 2 mesh2D -metis 2

# Serial cases
ElmerSolver Valid_CostRegSolver1.sif
# parallel cases
mpirun -np 2 ElmerSolver_mpi
```

