A few introductory runs to the double porous equivalent layer hydrological model

launching the $./MakeMSH command allow the build up of the mesh and the compilation of the code

  It is running as tha basal boundary condition of a small cube.

  The runs are : 

	SedOnly.sif : 	 Take into account only the sediment layer by using the "Dirichlet condition = FALSE" in the 
		    	 solver params. for the rest, the input is constant on all the domain, no flux on all but one
					 boundary where a water head is given. 
					 This allow an analytic solution : h(x)=(x^2-2Lx)*(-Q/2T) ; 
					 			L is the size of the domain, Q the input and the the sediment transmitivity

	Coupled.sif :    The same as the preceding but this time with an upper boundary for the sediment layer water
		    	 head. This settings allow to use the two layers of the system