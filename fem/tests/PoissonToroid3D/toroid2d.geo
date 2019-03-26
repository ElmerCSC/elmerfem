cl__1 = 1;
Point(1) = {2, 0, 0, 1};
Point(2) = {3, 0, 0, 1};
Point(3) = {2.5, 0.8, 0, 1};
Point(4) = {2.5, -0.8, 0, 1};
Spline(1) = {1, 4, 2, 3, 1};
Line Loop(3) = {1};
Plane Surface(3) = {3};
