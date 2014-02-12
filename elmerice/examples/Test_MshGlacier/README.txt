Compilation : elmerf90-nosh $ElmerIceLGGE/Meshers/MshGlacier.f90 $ELMER_HOME/lib/libelmersolver.dylib -o MshGlacier
ElmerGrid 1 2 name.grd
./MshGlacier
ElmerGrid 2 3 name 
