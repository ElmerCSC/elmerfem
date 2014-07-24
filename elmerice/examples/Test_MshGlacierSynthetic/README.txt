# Compilation
Extension of the library depends on system you are working on
# .so pour unix
----------------
elmerf90-nosh {path_to_elmerice}/Meshers/MshGlacierSynthetic.f90 fbed.f90 fsurf.f90 $ELMER_HOME/lib/libelmersolver.so -o MshGlacierSynthetic

# .dylib pour Mac
-----------------
elmerf90-nosh {path_to_elmerice}/Meshers/MshGlacierSynthetic.f90 fbed.f90 fsurf.f90 $ELMER_HOME/lib/libelmersolver.dylib -o MshGlacierSynthetic

# make the 1m thick mesh
ElmerGrid 1 2 mesh_A 
# Execution
./MshGlacierSynthetic
# For ElmerPost visualisation
ElmerGrid 2 3 mesh_A
