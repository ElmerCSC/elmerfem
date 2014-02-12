#!/bin/bash  

# Define the number of partitions (has to be a parallel run, i.e. np>1 !!)
export np=4 

# Test to run 

file=Adjoint_Beta.sif  # Control inverse method; optimisation of the slip coef.
#file=Adjoint_Beta_GradientValid.sif # Control inverse method; Compare Adjoint total derivative with finite difference
#file=Adjoint_Mu.sif    # Control inverse method; optimisation of the Viscosity
#file=Adjoint_Mu_GradientValid.sif # Control inverse method; Compare Adjoint total derivative with finite difference
#file=Robin_Beta.sif   # Robin inverse method; optimisation of the slip coef.
#file=Robin_Beta_GradientValid.sif # Robin inverse method; Compare Robin total derivative with finite difference

## make compilation and Mesh
make clean
make 
make Mesh

# Run the test
echo $file > ELMERSOLVER_STARTINFO
mpirun -n $np ElmerSolver_mpi
