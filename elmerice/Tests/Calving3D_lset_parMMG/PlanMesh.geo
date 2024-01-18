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
Physical Surface(105)= {106};

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
Mesh.CharacteristicLengthExtendFromBoundary = 0;//+

out[]=Extrude {0, 0, 1.0} {
  Surface{106}; Layers{5}; 
};

    /* surfaceVector contains in the following order:
    [0]	- front surface (opposed to source surface)
    [1] - extruded volume
    [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
    [3] - right surface (belonging to 2nd line in "Line Loop (6)")
    [4] - top surface (belonging to 3rd line in "Line Loop (6)")
    [5] - left surface (belonging to 4th line in "Line Loop (6)") */

n = #out[];
Printf("Extrude has returned %g elements", n);

Physical Volume(1) = {out[1]};
Physical Surface(106)= {out[0]};
Physical Surface(102)= {out[2]};
Physical Surface(103)= {out[3]};
Physical Surface(104)= {out[4]};
Physical Surface(101)= {out[5]};

// This mesh has boundaries as follows:
// 101 - calving front
// 102 - right sidewall
// 103 - inflow boundary
// 104 - left sidewall
// 105 - base
// 106 - surf