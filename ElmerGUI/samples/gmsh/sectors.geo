npts=12;
h=0.1;
r1=1.0;
r2=2.0;
h1=2*Pi*r1/(npts*10.0);
h2=2*Pi*r2/(npts*5.0);
For i In {1 : npts}
  theta = i * 2*Pi/npts;
  Point(i) = {r1 * Cos(theta), r1 * Sin(theta), 0.0, h1};
EndFor
For i In {1 : npts}
  theta = i * 2*Pi/npts;
  Point(i+npts) = {r2 * Cos(theta), r2 * Sin(theta), 0.0, h2};
EndFor
Point(2*npts+1) = {0.0, 0.0, 0.0, h1};
For i In {1 : npts}
  Line(i) = {i,i+npts};
EndFor
For i In {1 : npts-1}
  Circle(npts+i) = {i,2*npts+1,i+1};
EndFor
Circle(2*npts) = {npts,2*npts+1,1};
For i In {1 : npts-1}
  Circle(2*npts+i) = {npts+i,2*npts+1,npts+i+1};
EndFor
Circle(3*npts) = {2*npts,2*npts+1,npts+1};
For i In {1 : npts-1}
  Curve Loop(i) = {i, 2*npts+i, -(i+1), -(npts+i)};
  Plane Surface(i) = {i};
EndFor
//+
Curve Loop(npts) = {npts, 3*npts, -1, -(2*npts)};
Plane Surface(npts) = {npts};
//+
Curve Loop(npts+1) = {17, 18, 19, 20, 21, 22, 23, 24, 13, 14, 15, 16};
Plane Surface(npts+1) = {npts+1};

