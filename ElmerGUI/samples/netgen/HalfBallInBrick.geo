#
# a rectangle with half ball inside
#
algebraic3d

solid rectangle = orthobrick(0, 0, 0; 4, 2, 3);
solid ball = sphere (2.0, 0.0, 2.0; 0.5);

solid rest = rectangle and not ball;

# two sub-domains
tlo rest -col=[0,0,1];
