lc = 500.0;
Point(1)={1000.0,0.0,0.0,lc};
Point(2)={0.0,5000.0,0.0,lc};
Point(3)={5000.0,5000.0,0.0,lc};
Point(4)={4000.0,0.0,0.0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Physical Line(101) = {1};
Physical Line(102) = {2};
Physical Line(103) = {3};
Physical Line(104) = {4};

Line Loop(105) = {1,2,3,4};

Plane Surface(106) = {105};
Physical Surface(107) = {106};

Field[1] = Attractor;
Field[1].NNodesByEdge = 10;
Field[1].EdgesList = {4};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 100.0;
Field[2].LcMax = 500.0;
Field[2].DistMin = 250.0;
Field[2].DistMax = 1000.0;
Background Field = 2;
Mesh.CharacteristicLengthExtendFromBoundary = 0;