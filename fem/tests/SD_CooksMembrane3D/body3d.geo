Point(1) = {0, 0, 0, 2.0};
Point(2) = {48, 44, 0, 4.0};
Point(3) = {48, 60, 0, 2.0};
Point(4) = {0, 44, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Extrude {0,0,5} {
 Surface{6}; Layers{4};Recombine;
}
Physical Surface(29) = {27};
Physical Surface(30) = {19};
Physical Surface(31) = {28, 6};
Physical Surface(32) = {15, 23};
Physical Volume(33) = {1};
