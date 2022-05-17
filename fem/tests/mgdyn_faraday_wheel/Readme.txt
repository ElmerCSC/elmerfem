This is a special model of the Faraday wheel to test the mortar FE 
approximation. To this end the wheel is artificially broken into an inner part
and outer part which are meshed independently. Moreover, the vector potential
variable is made unique in an iterated manner (the nonlinear iteration is 
employed for this purpose) by performing the Helmholtz projection after 
solving the A-V solution candidate. For a similar verification case see the 
test ../mgdyn_faraday_disc/ and for a transient version the test
../mgdyn_faraday_disc_transient. The transient version of this problem with a
mortar BC is given in the file transient_with_projection.sif of the current 
directory.

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

The continuity of the A-V solution is enforced weakly over the artificial 
contact surface of the inner and outer part.
