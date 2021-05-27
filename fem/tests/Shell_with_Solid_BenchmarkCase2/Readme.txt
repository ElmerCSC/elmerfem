A verification case for assembling 3D solid and 2D shell equations
into a single linear system so that the coupled equations can 
be solved simultaneously.

This is a bending-dominated benchmark problem described in 
Pitkäranta et al. Shell deformation states and the finite element 
method: a benchmark study of cylindrical shells. Here the boundary
layer adjacent to the free edge is modelled in 3D and
the 2D shell model is used elsewhere. For a model using 
the shell equations only see the test 

      ../tests/Shell_BenchmarkCase2_Quad

The results of the coupled simulation have been verified
against using the pure shell equations.
