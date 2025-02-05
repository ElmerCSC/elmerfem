// Gmsh project created on Sat Dec 21 08:30:08 2024
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 20, 20, 20};
//+
Physical Surface("side", 13) = {6};
//+
Physical Surface("top", 14) = {4};
//+
Physical Surface("other", 15) = {1, 5, 2, 3};
//+
Physical Volume("body1", 16) = {1};
