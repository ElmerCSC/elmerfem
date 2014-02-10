% create teterousse.msh using gmsh
gmsh teterousse.geo -1 -2 

% convert teterousse.gmsh in an Elmer type mesh
ElmerGrid 14 2 teterousse.msh -autoclean

% Extrude vertically the mesh (1m thick)
% Get ExtrudeMesh from the svn 
% and compile it (cc ExtrudeMesh.c -o ExtrudeMesh -lm)
ExtrudeMesh teterousse WithOutCavity 14 1 1 0 0 0 0 

% Deform vertically using the surface and bedrock DEM
% Input data are in mesh_input.dat
elmerf90-nosh MshGlacier.f90 -o MshGlacier
./MshGlacierDEM

% Make a .ep to visualize in ElmerPost the mesh
ElmerGrid 2 3 WithOutCavity
