//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {1, 0, 0, 0.1};
//+
Point(3) = {0, 1, 0, 0.1};
//+
Point(4) = {-1, 0, 0, 0.1};
//+
Point(5) = {0, -1, 0, 0.1};
//+
Circle(1) = {2, 1, 3};
//+
Circle(2) = {3, 1, 4};
//+
Circle(3) = {4, 1, 5};
//+
Circle(4) = {5, 1, 2};
//+
Line(5) = {1, 2};
//+
Line(6) = {1, 3};
//+
Line(7) = {1, 4};
//+
Line(8) = {1, 5};
//+
Curve Loop(1) = {1, -6, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -7, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, -8, 7};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, -5, 8};
//+
Plane Surface(4) = {4};
