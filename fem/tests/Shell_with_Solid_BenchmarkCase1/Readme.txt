A verification case for assembling 3D solid and 2D shell equations
into a single linear system so that the coupled equations can 
be solved simultaneously.

This is a membrane-dominated benchmark problem described in 
Pitkäranta et al. Shell deformation states and the finite element 
method: a benchmark study of cylindrical shells. Here the boundary
layer adjacent to the clamped edge is modelled in 3D and
the 2D shell model is used elsewhere. For a model using 
the shell equations only see the test 

      ../tests/Shell_BenchmarkCase1_Quad

The coupling conditions are enforced quite near the clamped 
edge where the boundary layer is still strong in order to
test that the effects of bending and transverse shear deformation 
are handled correctly by the interface conditions.

The results of the coupled simulation have been verified
against using pure shell equations and pure 3D elasticity
equations.
