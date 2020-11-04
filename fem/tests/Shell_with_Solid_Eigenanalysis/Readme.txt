A verification case for assembling 3D solid and 2D shell equations
into a single linear system so that the eigenanalysis of the coupled 
equations can be performed.

This is a cylindrical shell problem with bending-dominated 
asymptotic behaviour. Here the boundary layer adjacent to the free 
edge is modelled in 3D and the 2D shell model is used elsewhere. 
The results of this simulation should be compared with those of 
the 2D shell model

       ../tests/Shell_BenchmarkCase2_Quad/eigenanalysis.sif

The eigenmodes computed with the different strategies
appear to be in good agreement. Note that only symmetric modes
are sought here.

For static analysis of the coupled equations see also the tests 

      ../tests/Shell_with_Solid_BenchmarkCase1/case.sif
      ../tests/Shell_with_Solid_BenchmarkCase2/case.sif
