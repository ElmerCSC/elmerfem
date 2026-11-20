// Unit cube [-1,1]^3 meshed with linear tetrahedra for the SplitFEMPoissonTets test.
// Create the mesh with:
//   gmsh -3 cube-tet.geo -o cube-tet.msh
SetFactory("OpenCASCADE");

Box(1) = {-1, -1, 0, 2, 2, 0.2};

// OpenCASCADE numbers the box faces x-min, x-max, y-min, y-max, z-min, z-max.
Physical Surface("xmin", 1) = {1};
Physical Surface("xmax", 2) = {2};
Physical Surface("ymin", 3) = {3};
Physical Surface("ymax", 4) = {4};
Physical Surface("zmin", 5) = {5};
Physical Surface("zmax", 6) = {6};
Physical Volume("cube", 1) = {1};

Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.05;
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 2.2;
