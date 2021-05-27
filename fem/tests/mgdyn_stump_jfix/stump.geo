Point(1) = {0, 0, 0, 0.02};
Point(2) = {0.1, 0, 0, 0.02};
Point(3) = {0, 0.1, 0, 0.02};
Point(4) = {-0.1, 0, 0, 0.02};
Point(5) = {0, -0.1, 0, 0.02};
Circle(1) = {3, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 3};
Line Loop(5) = {4, 1, 2, 3};
Ruled Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6};
}
Point(19) = {-1, -1, -1, 0.2};
Point(20) = {1, -1, -1, 0.2};
Point(21) = {1, 1, -1, 0.2};
Point(22) = {-1, 1, -1, 0.2};
Line(29) = {22, 21};
Line(30) = {21, 20};
Line(31) = {19, 19};
Line(32) = {20, 19};
Line(33) = {22, 19};
Extrude {0, 0, 3} {
  Line{33, 29, 30, 32};
}
Line Loop(50) = {34, -46, -42, -38};
Plane Surface(51) = {50};
Line Loop(52) = {33, -32, -30, -29};
Plane Surface(53) = {52};
Surface Loop(54) = {37, 53, 49, 51, 45, 41};
Surface Loop(55) = {28, 15, 6, 19, 23, 27};
Volume(56) = {54, 55};
