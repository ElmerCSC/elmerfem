c = 1.0;
h1 = c*0.1;
h2 = c*0.3;
r1 = 1.0;
r2 = 1.2;
r3 = 3.0;
r4 = r2;
h3 = c*h1;  //(h1+h2)/2;
Point(1) = {0, 0, 0, h1};
Point(2) = {r1, 0, 0, h1};
Point(3) = {0, r1, 0, h1};
Point(4) = {-r1, 0, 0, h1};
Point(5) = {0, -r1, 0, h1};
Point(6) = {r2, 0, 0, h3};
Point(7) = {0, r2, 0, h3};
Point(8) = {-r2, 0, 0, h3};
Point(9) = {0, -r2, 0, h3};

Point(10) = {r3, 0, 0, h2};
Point(11) = {0, r3, 0, h2};
Point(12) = {-r3, 0, 0, h2};
Point(13) = {0, -r3, 0, h2};

Point(14) = {r4, 0, 0, h3};
Point(15) = {0, r4, 0, h3};
Point(16) = {-r4, 0, 0, h3};
Point(17) = {0, -r4, 0, h3};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
//+
Line(9) = {4, 2};
//+
Curve Loop(1) = {2, 9, 1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, -9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, 7, 8, 5};
//+
Curve Loop(4) = {2, 3, 4, 1};
//+
Plane Surface(3) = {3, 4};
//+
Circle(10) = {11, 1, 10};
//+
Circle(11) = {10, 1, 13};
//+
Circle(12) = {13, 1, 12};
//+
Circle(13) = {12, 1, 11};
//+
Circle(14) = {15, 1, 14};
//+
Circle(15) = {14, 1, 17};
//+
Circle(16) = {17, 1, 16};
//+
Circle(17) = {16, 1, 15};
//+
Curve Loop(5) = {13, 10, 11, 12};
//+
Curve Loop(6) = {17, 14, 15, 16};
//+
Plane Surface(4) = {6, 5};
