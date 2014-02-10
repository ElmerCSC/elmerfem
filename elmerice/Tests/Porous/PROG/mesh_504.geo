lc = 0.10 ;
Nz = 10 ; 
Point(1) = {0.0,0.0,0.,lc};
Point(2) = {1.0,0.0,0.,lc};
Point(3) = {1.0,1.0,0.,lc};
Point(4) = {0.0,1.0,0.,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {1:4};
Plane Surface(6) = {5};


out[] = Extrude {0,0,1} {
//  Line{1}; Layers{ {Nz},{1} } ;  Recombine ;
  Line{1}; Layers{ {Nz},{1} } ;  
} ;
Physical Surface(1) = {out[1]};


out[] = Extrude {0,0,1} {
//  Line{2}; Layers{ {Nz},{1} } ; Recombine ;
  Line{2}; Layers{ {Nz},{1} } ; 
} ;
Physical Surface(2) = {out[1]};

out[] = Extrude {0,0,1} {
//  Line{3}; Layers{ {Nz},{1} } ; Recombine ;
  Line{3}; Layers{ {Nz},{1} } ; 
} ;
Physical Surface(3) = {out[1]};

out[] = Extrude {0,0,1} {
//  Line{4}; Layers{ {Nz},{1} } ; Recombine ;
  Line{4}; Layers{ {Nz},{1} } ; 
} ;
Physical Surface(4) = {out[1]};


Physical Surface(5) = {6};

out[] =  Extrude {0,0,1} {
//  Surface{6}; Layers{ {Nz},{1} } ; Recombine ; 
  Surface{6}; Layers{ {Nz},{1} } ;  
} ;  

Physical Surface(6) = {out[0]};

vol = out[1] ;

//out[] = Extrude {0,0,100} {
//  Point{1}; Layers{ {Nz},{1} } ;  
//} ;
//Physical Point(7) = {out[0]} ;

Physical Volume(8) = {vol} ;

