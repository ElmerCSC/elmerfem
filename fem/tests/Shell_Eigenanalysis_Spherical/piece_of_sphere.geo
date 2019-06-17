//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {1, 0, 0, 0.1};
//+
Point(3) = {0.2, 0.8, 0, 0.01};
//+
Point(4) = {0, 0.8, 0, 0.01};
//+
Point(5) = {0, 1, 0, 0.01};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 5};
//+
Line(4) = {5, 4};
//+
Line(5) = {4, 1};
//+
Line(6) = {3, 4};
//+
Curve Loop(1) = {1, 2, 6, 5};
Curve Loop(2) = {6, -4, -3};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Transfinite Line{1} = 41;
Transfinite Line{6} = 41;
Transfinite Line{2} = 33;
Transfinite Line{5} = 33;
Transfinite Surface {1} = {1,2,3,4};
Recombine Surface{1};
//+
Physical Curve("equator") = {1};
//+
Physical Curve("meridian_1") = {2, 3};
//+
Physical Curve("meridian_2") = {4, 5};
//+
Physical Surface("mid_surface") = {1, 2};
