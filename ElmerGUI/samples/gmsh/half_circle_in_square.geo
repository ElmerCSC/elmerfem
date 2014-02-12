// characteristic mesh parameters
c1 = 0.1;
c2 = 0.02;
// points for the rectangle
Point(1) = {0, -1, 0, c1};
Point(2) = {0, 1, 0, c1};
Point(3) = {1, 1, 0, c1};
Point(4) = {1, -1, 0, c1};
// points for the half circle
Point(5) = {0, 0, 0, c2};
Point(6) = {0, 0.1, 0, c2};
Point(7) = {0.1, 0, 0, c2};
Point(8) = {0, -0.1, 0, c2};
// define the lines
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 8};
Line(5) = {6, 2};
Circle(6) = {8, 5, 7};
Circle(7) = {7, 5, 6};
Line(8) = {8, 5};
Line(9) = {5, 6};
// define the surfaces
Line Loop(10) = {5, 1, 2, 3, 4, 6, 7};
Plane Surface(11) = {10};
Line Loop(12) = {6, 7, -9, -8};
Plane Surface(13) = {12};
