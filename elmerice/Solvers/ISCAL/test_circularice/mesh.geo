// Mesh construction:
//------------------
//Create a mesh in gmsh using the file mesh.geo : q and dphi determines mesh resolution.
//Make the mesh into elmerformat by using ElmerGrid 14 2 -autoclean mesh.msh



cl1 = 100000;

pi=3.14;

dphi=20;
n=dphi*2-1; //number of points circular direction
q=12; //number of points in radial direction

Point(1)={0, 0, 0, cl1};

rinner=450000;
router=750000;
dr=(router-rinner)/(q-1);

r[0]=rinner;
For k In {1:q-1:1}
r[k]=r[k-1]+dr;
EndFor
r[q]=router;

phi=0;

//Coordinates and points
For k In {0:q-1:1}
x[k]=r[k]*Cos(phi);
y[k]=r[k]*Sin(phi);
Point(2+k*(n+1)) = {x[k], y[k], 0, cl1};
EndFor

//Radial line
i=0;
For k In {0:q-2:1}
Line(10000+k*(n+1)+i) = {2+k*(n+1)+i,2+(k+1)*(n+1)+i};
EndFor

v=0;

For i In {1:n:1}

//Angles
phi=i*Pi/dphi;

//Coordinates and points
For k In {0:q-1:1}
x[k]=r[k]*Cos(phi);
y[k]=r[k]*Sin(phi);
Point(2+k*(n+1)+i) = {x[k], y[k], 0, cl1};
EndFor

//Circle arches
For k In {0:q-1:1}
	Circle(k*(n+1)+i)={2+k*(n+1)+(i-1),1,2+k*(n+1)+i};
	If(k==q-1)
		thelines[v]=k*(n+1)+i;
	v=v+1;
	EndIf
EndFor

//Radial lines
For k In {0:q-2:1}
Line(10000+k*(n+1)+i) = {2+k*(n+1)+i,2+(k+1)*(n+1)+i};
EndFor

EndFor


//Plane Surfaces, made out of: {inner arch, left line, outer arch, right line}

l=0;
For i In {0:n-1:1}
For k In {0:q-2:1}

Line Loop(14000+k*(n+1)+i) = {1+(k*(n+1)+i), 1+(10000+k*(n+1)+i), -(1+(k+1)*(n+1)+i), -(10000+k*(n+1)+i)};
Plane Surface(16000+k*(n+1)+i) = {14000+k*(n+1)+i};
Recombine Surface {16000+k*(n+1)+i};
EndFor
EndFor

For i In {0:n-1:1}
For k In {0:q-2:1}
Transfinite Surface{16000+k*(n+1)+i}={2+k*(n+1)+(i+1),2+k*(n+1)+i,2+(k+1)*(n+1)+(i+1),2+(k+1)*(n+1)+i};

thesquares[l]=16000+k*(n+1)+i;
l=l+1;

EndFor
EndFor


//--------------------------
//Closing the circle. 
//--------------------------

//Arches made of right point (i=n), center, left point (i=0)

//Create arches (i=n+1):
i=n+1;
For k In {0:q-1:1}
Circle(70000+k)={2+k*(n+1)+n,1,2+k*(n+1)+0};
If (k==q-1)
thelines[v]=k*(n+1)+i;
v=v+1;
EndIf
EndFor



//Plane surfaces (i=n) made out of {inner arch (new), left line (i=0), outer arch (new), right line (i=n) }

For k In {0:q-2:1}

Line Loop(80000+k) = {70000+k,10000+k*(n+1)+0, -(70001+k), -(10000+k*(n+1)+n)};
Plane Surface(90000+k) = {80000+k};
Recombine Surface {90000+k};

Transfinite Surface{90000+k}={2+k*(n+1)+0,2+k*(n+1)+n,2+(k+1)*(n+1)+0,2+(k+1)*(n+1)+n};

thesquares[l]=90000+k;
l=l+1;

EndFor


//inner circle plane surface

hh[]={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 
41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 
61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 
82,83,84,85,86,87,88,89,90,
91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,
111,112,113,114,115,116,117,118,119,
120,121,122,123,124,125,126,127,128,129};

Line Loop(90000+q-1)={hh[0]:hh[n-1],70000};

Plane Surface(90000+q) = {90000+q-1};

//Combine to one physical surface with one physical boundary

Physical Line(110000) = {thelines[]};

Physical Surface(150000) = {thesquares[],90000+q};

