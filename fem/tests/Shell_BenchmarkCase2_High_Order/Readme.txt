Solve nonlinear shell equations in the case of a cylindrical benchmark 
problem with bending-dominated asymptotic behaviour in the linear regime. 
The linearized problem has been described in Pitkäranta et al. Shell 
deformation states and the finite element method: a benchmark study of 
cylindrical shells. Computer Methods in Applied Mechanics and Engineering 
1995. 128:81-121. As opposed to the reference, here the shell
problem is treated as nonlinear, but the shell is subject to a small load 
so that nonlinear effects are not yet significant.

Here the problem is solved by using a high-order discretization over a 2-D
domain. This option needs a special command "Skip Surface Reconstruction = 
True" so that instead of using a physical surface mesh the surface 
parametrization is described by the subroutine SurfaceBasis within the solver
code ShellSolver.F90. Working without a physical surface mesh is limited to
some special geometries which allow a parametrization by lines of curvature
coordinates. There are other tests for general shell modelling based on
physical surface meshes.
