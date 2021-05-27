The solution of a magnetostatic problem by using the consistently regularized 
formulation (to ensure a unique vector potential solution). See also the 
alternate test mgdyn_bh_gauge2 for using a more effective linear solver.  

The consistently regularized formulation employs the scalar variable as 
an augmented Lagrangian variable to impose the divergence-free constraint on A. 
The uniqueness of the vector potential can also be enforced in an iterated
manner without introducing a coupled system to solve A and the Lagrange 
multiplier. For this alternative see the file with_projection.sif in this 
directory.
