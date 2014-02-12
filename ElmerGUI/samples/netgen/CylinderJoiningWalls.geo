#
# a rectangle with half ball inside
#
algebraic3d

solid rect1 = orthobrick(-3.0, -3.0, -3.0; -2.0, 3, 3);
solid rect2 = orthobrick(2.0, -3, -3; 3, 3, 3);
solid rect3 = orthobrick(-3, -3, -3; 3, 3, 3);
solid cyl = cylinder(-3.0, 0.0, 0.0; 3.0, 0.0, 0.0; 1.0);

solid rest = rect1 or rect2 or (cyl and rect3);

# two sub-domains
tlo rest -col=[0,0,1];
