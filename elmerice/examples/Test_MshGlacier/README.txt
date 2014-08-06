# Compilation in mismip or Helheim directory #
----------------------------------------------

Extension of the library depends on system you are working on

# .so pour unix
----------------
elmerf90-nosh {path_to_elmerice}/Meshers/MshGlacier.f90 $ELMER_HOME/lib/libelmersolver.so -o MshGlacier

# .dylib pour Mac
-----------------
elmerf90-nosh {path_to_elmerice}/Meshers/MshGlacier.f90 $ELMER_HOME/lib/libelmersolver.dylib -o MshGlacier


ElmerGrid 1 2 name.grd
./MshGlacier

# Make a .ep to visualize in ElmerPost the mesh
ElmerGrid 2 3 name

