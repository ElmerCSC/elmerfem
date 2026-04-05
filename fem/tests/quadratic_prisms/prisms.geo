Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Set mesh options for second order prisms
Mesh.SecondOrderIncomplete = 0;

// Extrude the surface to create a volume
Extrude {0, 0, 1} {
  Surface{1}; Layers {5}; Recombine;
}
