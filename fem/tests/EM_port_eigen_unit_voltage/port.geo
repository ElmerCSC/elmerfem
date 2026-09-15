//+
Point(1) = {0, -25, 0, 2.0};
//+
Point(2) = {50, -25, 0, 2.0};
//+
Point(3) = {50, -10.2, 0, 2.0};
//+
Point(4) = {50, 25, 0, 2.0};
Point(5) = {0, 25, 0, 2.0};
Point(6) = {0, -10.2, 0, 2.0};

Point(7) = {14, -10.2, 0, 0.05};
Point(8) = {36, -10.2, 0, 0.05};
Point(9) = {36, -10.0, 0, 0.05};
Point(10) = {14, -10.0, 0, 0.05};




//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 7};
Line(11) = {6, 7};
Line(12) = {8, 3};

//+
Curve Loop(1) = {1, 2, -12, -7, -11, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, 4, 5, 11, -10, -9, -8, 12};
//+
Plane Surface(2) = {2};
//+
Physical Surface("substrate", 13) = {1};
//+
Physical Surface("air", 14) = {2};
//+
Physical Curve("outersurf", 15) = {1, 2, 3, 4, 5, 6};
//+
Physical Curve("signal", 16) = {7, 8, 9, 10};

