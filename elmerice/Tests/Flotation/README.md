# Test Sub-Element GL parameterization for Flotation

## Description
Simple ice shelf with uniform thickness.

Test Sub-Element GL parameterizations to compute the ice fractions
- Case 1 = NoSEP : using standard number of IP
- Case 2 = SEP2 : Adaptive Integration with Split FEM
- Case 3 = SEP3 : increasing number of IP

Analytically the GL is located at x=170km, the grounded area is then 170.0e3x50.0e3=8.5e9 m^2 and is obtained exactly with SEP2 

## Running the test
1. Generate the mesh :
	> ElmerGrid 1 2 mesh.grd
2. Run the test:
	> ElmerSolver 

 ## Results established:

- 08.04.2026
- Fabien Gillet-Chaulet, IGE
