This is a simple verification case to check that a Faraday disc creates
the expected potential difference between the rim and and the axis.
As explained in the following, this case is also instructive to understand 
how Elmer treats Maxwell's equations in the case of a moving body.

Here a metallic disc rotates and is subject to a uniform magnetic field 
parallel to the disc axis. When the disc is put to rotate, the mobile charges
of the conductor start to feel a radial force due to the Lorentz term v x B
and hence accumulate on the insulated rim surface. As a consequence of charge
separation a radial electric field E is produced. The process of charge 
separation continues as long as steady state is reached such that the force on
the charges which is produced by the Lorentz term v x B balances the force 
which is produced by E. That is, in the steady state we have 

              E'(x) = E(x) + v(x) x B(x) = 0, 

with E'(x) being the spatial field that describes the force experienced by 
a charge (see Elmer Models Manual for Maxwell's equations in terms of E'). 
Note that the solution is expressed in the spatial coordinates and attention 
is not fixed on particular material points (this is known as the Euler 
description in the continuum mechanics).  

The potential difference between the centre and the rim is known to be 
1/2 w a^2 |B|, where a is the radius of the disc and w is the angular 
velocity. Elmer generates E'(x) as the electric field variable and we expect 
to obtain E' = 0 up to the discretization error.
