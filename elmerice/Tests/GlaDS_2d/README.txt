To execute the test:
--------------------

#################
# Run on 1 proc
#################


# Run 
mpirun -np 1 ElmerSolver_mpi glads_2d.sif 


#################
# Run on 2 procs
#################

# Create 2  partitions
ElmerGrid  2 2 mesh_B5  -partition  2 1 1

# Distribute Moulins
python makemoulin.py --meshdir mesh_B5  --moulin B5_M.xy  --partition 2

# Run
mpirun -np 2 ElmerSolver_mpi glads_2d.sif 



###############################
# Notes: Create Mesh from .geo
###############################


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
Revision 926851c7 


