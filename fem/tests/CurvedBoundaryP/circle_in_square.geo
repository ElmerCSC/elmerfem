// characteristic mesh parameters
h = 1.0;
r = 0.5;
c1 = 0.5;
c2 = 0.5;
// points for the rectangle
Point(1) = {-h, -h, 0, c1};
Point(2) = {-h, h, 0, c1};
Point(3) = {h, h, 0, c1};
Point(4) = {h, -h, 0, c1};
// points for the half circle
Point(5) = {0, 0, 0, c2};
Point(6) = {0, r, 0, c2};
Point(7) = {r, 0, 0, c2};
Point(8) = {0, -r, 0, c2};
Point(9) = {-r, 0, 0, c2};
// define the lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Circle(5) = {7, 5, 6};
//+
Circle(6) = {6, 5, 9};
//+
Circle(7) = {9, 5, 8};
//+
Circle(8) = {8, 5, 7};
//+
Curve Loop(1) = {6, 7, 8, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {1, 2, 3, 4};
//+
Plane Surface(2) = {1, 2};
