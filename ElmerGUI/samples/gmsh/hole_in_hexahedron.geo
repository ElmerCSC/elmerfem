// the square
Point(1) = {0, 0, 0, 0.1};
Point(2) = {0, 1, 0, 0.1};
Point(3) = {1, 1, 0, 0.1};
Point(4) = {1, 0, 0, 0.1};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
// the circle in the centre
Point(5) = {0.5, 0.5, 0, 0.1};
Point(6) = {0.5, 0.7, 0, 0.1};
Circle(5) = {6, 5, 6};
// creating the surfaces
Line Loop(6) = {2, 3, 4, 1};
Line Loop(7) = {5};
Plane Surface(8) = {6, 7};
Plane Surface(9) = {7};
//extruding them
Extrude {0, 0, 1} {
  Surface{9, 8}; Layers{10};
}
// bodies and boundaries
Physical Volume(44) = {1};
Physical Volume(45) = {2};
Physical Surface(46) = {9};
Physical Surface(47) = {16};
Physical Surface(48) = {8};
Physical Surface(49) = {43};
Physical Surface(50) = {15};
Physical Surface(51) = {34};
Physical Surface(52) = {26};
Physical Surface(53) = {30};
Physical Surface(54) = {38};
