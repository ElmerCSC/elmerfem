A benchmark problem where an open hemispherical shell is subject to point
loads. The reference results can be found in Sze KY, Liu XH, Lo SH. Popular 
benchmark problems for geometric nonlinear analysis of shells. Finite 
Elements in Analysis and Design 2004, 40(11):1551-1569.

Here the problem is solved by using a high-order discretization over a 2-D
domain. This option needs a special command "Skip Surface Reconstruction = 
True" so that instead of using a physical surface mesh the surface 
parametrization is described by the subroutine SurfaceBasis within the solver
code ShellSolver.F90. Working without a physical surface mesh is limited to
some special geometries which allow a parametrization by lines of curvature
coordinates. There are other tests for general shell modelling based on
physical surface meshes.
