//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {1, 0, 0, 0.1};
//+
Point(3) = {1, 1, 0, 0.1};
//+
Point(4) = {0, 1, 0, 0.1};
//+
Point(5) = {0, 0, 1, 0.1};
//+
Point(6) = {1, 0, 1, 0.1};
//+
Point(7) = {1, 1, 1, 0.1};
//+
Point(8) = {0, 1, 1, 0.1};
//+
Point(9) = {0.5, 0.25, 0, 0.1};
//+
Point(10) = {0.5, 0.75, 0, 0.1};
//+
Point(11) = {0.5, 0.25, 1, 0.1};
//+
Point(12) = {0.5, 0.75, 1, 0.1};
//+
Line(1) = {5, 6};
//+
Line(2) = {6, 7};
//+
Line(3) = {7, 8};
//+
Line(4) = {8, 5};
//+
Line(5) = {6, 2};
//+
Line(6) = {7, 3};
//+
Line(7) = {12, 11};
//+
Line(8) = {11, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 12};
//+
Line(11) = {1, 4};
//+
Line(12) = {4, 8};
//+
Line(13) = {5, 1};
//+
Line(14) = {1, 2};
//+
Line(15) = {2, 3};
//+
Line(16) = {3, 4};
//+
Curve Loop(1) = {14, -5, -1, 13};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {15, -6, -2, 5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, -12, -16, -6};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 12, 4, 13};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 1, 2, 3};
//+
Curve Loop(6) = {7, 8, 9, 10};
//+
Line(17) = {5, 11};
//+
Line(18) = {11, 6};
//+
Line(19) = {12, 8};
//+
Line(20) = {12, 7};
//+
Line(21) = {10, 4};
//+
Line(22) = {10, 3};
//+
Line(23) = {9, 2};
//+
Line(24) = {9, 1};
//+
Curve Loop(7) = {1, -18, -17};
//+
Plane Surface(5) = {7};
//+
Curve Loop(8) = {18, 2, -20, 7};
//+
Plane Surface(6) = {8};
//+
Curve Loop(9) = {20, 3, -19};
//+
Plane Surface(7) = {9};
//+
Curve Loop(10) = {19, 4, 17, -7};
//+
Plane Surface(8) = {10};
//+
Curve Loop(11) = {16, -21, 22};
//+
Plane Surface(9) = {11};
//+
Curve Loop(12) = {21, -11, -24, 9};
//+
Plane Surface(10) = {12};
//+
Curve Loop(13) = {24, 14, -23};
//+
Plane Surface(11) = {13};
//+
Curve Loop(14) = {23, 15, -22, -9};
//+
Plane Surface(12) = {14};
//+
Plane Surface(13) = {6};
//+
Curve Loop(15) = {17, 8, 24, -13};
//+
Plane Surface(14) = {15};
//+
Curve Loop(16) = {18, 5, -23, -8};
//+
Plane Surface(15) = {16};
//+
Curve Loop(17) = {20, 6, -22, 10};
//+
Plane Surface(16) = {17};
//+
Curve Loop(18) = {19, -12, -21, 10};
//+
Plane Surface(17) = {18};
//+
Surface Loop(1) = {3, 7, 9, 17, 16};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {1, 11, 5, 15, 14};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {2, 12, 6, 15, 13, 16};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {14, 13, 17, 4, 10, 8};
//+
Volume(4) = {4};
