Test Case for transient mesh adaptation.

Use the ThicknessSolver to solve the problem of bodies in rotation under a divergence free velocity field.
This is a common test case to test codes solving a pure transport equation.


The mesh adaptation is performed using the algorithm described in:
	Alauzet, F., Frey, P.J., George, P.L., Mohammadi, B., 2007. 3D transient fixed point mesh adaptation for time-dependent problems: Application to CFD simulations. Journal of Computational Physics 222, 592–623. https://doi.org/10.1016/j.jcp.2006.08.012

The algorithm is implmented as a bash script avaialble in $SOURCE_DIR/elmerice/Solvers/MeshAdaptation_2D/Script_Transient.sh
Documentation can be found @: http://elmerice.elmerfem.org/wiki/doku.php?id=mesh:meshadaptation

F. Gillet-Chaulet
IGE, Grenoble
Nov. 2017
