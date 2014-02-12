# Make the mesh 
elmerf90-nosh MshGlacierSynthetic.f90 fbed.f90 fsurf.f90 -o MshGlacierSynthetic
ElmerGrid 1 2 mesh 
./MshGlacierSynthetic
