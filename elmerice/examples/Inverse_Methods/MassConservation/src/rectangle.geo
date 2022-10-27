// Generic gmsh .geo file to mesh a rectangular domain
// parameters: resolution and rectangle length and width
//  (can be changed in line using e.g. "-setnumber lc [VALUE]")
DefineConstant[ lc = {1.0e3, Name "resolution"},
                x = {200.0e3, Name "length"},
		y = {50.0e3, Name "width"}];

////////////////////////////////////////////////////////
Mesh.Algorithm=5;
// Points
Point(1) = {0, 0, 0, lc};
Point(2) = {x, 0, 0, lc};
Point(3) = {x, y, 0, lc};
Point(4) = {0, y, 0, lc};
// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Physical 
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
//+
Physical Surface(5) = {1};
