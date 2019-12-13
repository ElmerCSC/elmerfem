This is an axially symmetric test case where the UMAT subroutine is called to 
obtain the definition of the constitutive law of St Venant-Kirchhoff type. 
The case is derived from the test .../tests/TankAst (see this root case for the
origin of the test problem), but here the problem is handled as nonlinear and 
using the local UMAT subroutine to define the constitutive law. Since the load
is not severe, the results are essentially the same as obtained by using the 
linear model of this problem (and without UMAT). 

The author of this version: Mika Malinen
