This is a simple 2-D test case where the UMAT subroutine is called to obtain 
the definition of the basic isotropic constitutive law (the case of plane 
strain). The small strain tensor is employed, so the nonlinear iteration should
terminate after the first solution step. To increment the load without inertia
effects, the simulation type is here selected to be "scanning". The case has 
been checked to yield the same result as obtained by using the linear elasticity
solver (StressSolve).
