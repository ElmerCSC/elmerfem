ElmerGrid 14 2 mesh.msh -autoclean

elmerf90 -o linkFS linkFS.f90

ElmerSolver iscal.sif
