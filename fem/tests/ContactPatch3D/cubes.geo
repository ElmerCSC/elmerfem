Point(1) = {0, 0, 0, 0.5};

N = 4;
N2 = 3;

Extrude {1, 0, 0} {
  Point{1};
  Layers{N};
  Recombine;
}
Extrude {0, 1, 0} {
  Line{1};
  Layers{N};
  Recombine;
}

Extrude {0, 0, 1} {
  Surface{5};
  Layers{N};
  Recombine;
}

Point(15) = {0, 0, 1.00001, 0.5};

Extrude {1, 0, 0} {
  Point{15};
  Layers{N2};
  Recombine;
}
Extrude {0, 1, 0} {
  Line{28};
  Layers{N2};
  Recombine;
}
Extrude {0, 0, 1} {
  Surface{32};
  Layers{N2};
  Recombine;
}

Physical Volume(55) = {1}; // lower block
Physical Volume(56) = {2}; // upper block
Physical Surface(57) = {5}; // bottom of lower block
Physical Surface(58) = {27}; // top of lower block
Physical Surface(59) = {32}; // bottom of upper block
Physical Surface(60) = {54}; // top of upper block
Physical Surface(61) = {18, 45}; // symmetry y-z
Physical Surface(62) = {49, 22}; // symmetry x-z
