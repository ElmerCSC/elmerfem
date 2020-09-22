To execute the test:
--------------------

#################
# Run on 1 proc
#################

# Run 
mpirun -np 1 ElmerSolver_mpi glads_3dint.sif 


#################
# Run on 2 procs
#################

# Replace UMFPACK with MUMPS in the sif
# file if you are using  multiple partitions

# Mesh with 2 partitions (mesh_B5/partitioning.2)
# including moulins nodes is provided with this test 


# Create 2  partitions
ElmerGrid  2 2 mesh_B5  -partition  2 1 1

# Distribute Moulins
python makemoulin.py --meshdir mesh_B5  --moulin B5_M.xy  --partition 2

# Run
mpirun -np 2 ElmerSolver_mpi glads_3dint.sif 


###############################
# Notes: Create Mesh from .geo
###############################

# Need gmsh version 3.0.6 or less
# Otherwise you will not be able to create the same mesh 

gmsh -1 -2 mesh_B5.geo
ElmerGrid 14 2 mesh_B5.msh -autoclean

# Then distribute Moulin Nodes among the partitions
python makemoulin.py --meshdir mesh_B5  --moulin B5_M.xy  --partition NUMBER_OF_PARTITIONS



Results established:
------------------
22.12.2019
Mondher CHEKKI  IGE
Gricad Cluster: Dahu
Compiler version:  intel/intelmpi version 18.0.5
                   gcc/gfortran version 6.3.0
Revision bddcb7f2


