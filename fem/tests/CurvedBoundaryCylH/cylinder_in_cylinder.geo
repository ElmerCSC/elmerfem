// characteristic mesh parameters
a = 0.5;
b = 1.0;
c1 = 0.5;
c2 = 0.5;
// points for the outer circle
Point(1) = {0, b, 0, c1};
Point(2) = {b, 0, 0, c1};
Point(3) = {0, -b, 0, c1};
Point(4) = {-b, 0, 0, c1};
// points for the inner circle
Point(5) = {0, 0, 0, c2};
Point(6) = {0, a, 0, c2};
Point(7) = {a, 0, 0, c2};
Point(8) = {0, -a, 0, c2};
Point(9) = {-a, 0, 0, c2};
// define the lines
Circle(1) = {1, 5, 2};
Circle(2) = {2, 5, 3};
Circle(3) = {3, 5, 4};
Circle(4) = {4, 5, 1};
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
//+
Extrude {0, 0, 0.5} {
  Surface{1}; Surface{2}; 
}
