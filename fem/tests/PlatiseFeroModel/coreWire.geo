SetFactory("OpenCASCADE");

Cylinder(1) = {0, 0, -0.15, 0, 0, 0.3,  0.100, 2*Pi};
Cylinder(2) = {0, 0, -0.15, 0, 0, 0.3,  0.002, 2*Pi};
Torus(3)    = {0, 0, 0, 0.02, 0.003, Pi};

Rotate {{0, 0, 1}, {0, 0, 0}, Pi} {
  Duplicata { Volume{3}; }
}

BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4}; Delete; }

Physical Surface("Sky", 16) = {9, 7, 8};
Physical Surface("Wire End 1", 17) = {5};
Physical Surface("Wire End 2", 18) = {6};
Physical Surface("Wire Surface", 19) = {4};
Physical Surface("Core Surface", 20) = {11, 10};
Physical Surface("Core Flux 1", 21) = {12};
Physical Surface("Core Flux 2", 22) = {13};

Physical Volume("Air", 23) = {5};
Physical Volume("Wire", 24) = {2};
Physical Volume("Core", 25) = {3, 4};
