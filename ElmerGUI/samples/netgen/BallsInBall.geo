#
# a large ball with four circular holes 
#
algebraic3d

solid smallballs = sphere (-0.4, -0.3, -0.2; 0.26)
           or sphere (0.4, -0.3, -0.2; 0.18)
           or sphere (0.0, 0.5, -0.2; 0.17)
           or sphere (0.0, 0.0, 0.4; 0.23);

solid bigball = sphere (0.0, 0.0, 0.0; 1.0);

solid rest = bigball and not smallballs;

# two sub-domains
tlo smallballs -col=[1,0,0];
tlo rest -col=[0,0,1] -transparent;
