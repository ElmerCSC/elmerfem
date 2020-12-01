# Solvers

This set of documentation presents the various solvers developed to solve glaciological applications using the finite element code Elmer. A solver is a Fortran subroutine that solves a differential equation on a given domain for all the nodal variables belonging to this domain. Very often, this results in solving a linear system. In certain cases, when the equation is only local, this can be done more simply by updating and iterating over all the nodes of the domain.

Solvers in the Elmer/Ice package can be called using the object file name ElmerIceSolvers:
`Procedure =  File "ElmerIceSolvers" "NameSolver"`

where NameSolver is the name of the solver you want to execute. The source code of the solvers of the Elmer/Ice package can be obtained from the Elmer svn in /trunk/elmerice/Solvers/.

Non-distributed solvers must be linked to libelmersolver.so, for which a wrapper script elmerf90 exists:
`elmerf90 MySolver.f90 -o MySolver`

In the SIF file (Solver section), the solver has to point on this object file MySolver:
`Procedure =  File "MySolver" "NameProcedureInMySolver"`

where NameProcedureInMySolver is the name of the fortran procedure in the file MySolver.f90.
