Solve nonlinear shell equations of Reissner-Naghdi type when the shell 
solver reads the director data as an ordinary Elmer solver variable. This
case can also be run in parallel, while parallel execution is not currently
possible if the data file format mesh.director or mesh.elements.data is used
to specify the director (currently there is no parallel file formats for these
files). 

This case employs a slightly modified version of the solver NormalSolver
(the source code is contained in CylinderNormalSolver.F90).
The original mesh files located in the mesh directory give nodes which are on
a cylindrical surface. The known relation between the global coordinates and
the director vector is used to solve the Elmer variable 'Director' that
represents the shell director data.

A cylindrical benchmark problem with membrane-dominated asymptotic behaviour 
in the linear regime is solved. The linearized problem has been described in 
Pitkäranta et al. Shell deformation states and the finite element method: 
a benchmark study of cylindrical shells. Computer Methods in Applied Mechanics
and Engineering 1995. 128:81-121. Here the shell is subject to a moderate 
load so that nonlinear effects start to be significant. 
