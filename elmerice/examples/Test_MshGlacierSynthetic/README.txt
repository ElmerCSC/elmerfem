# Compilation
elmerf90-nosh MshGlacierSynthetic.f90 fbed.f90 fsurf.f90 -o MshGlacierSynthetic
# make the 1m thick mesh
ElmerGrid 1 2 mesh_A 
# Execution
./MshGlacierSynthetic
# For ElmerPost visualisation
ElmerGrid 2 3 mesh_A
