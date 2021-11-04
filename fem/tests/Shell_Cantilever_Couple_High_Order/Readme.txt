A benchmark problem where a cantilever is subject to a moment load at an end.
The reference results can be found in Sze KY, Liu XH, Lo SH. Popular 
benchmark problems for geometric nonlinear analysis of shells. Finite 
Elements in Analysis and Design 2004, 40(11):1551-1569.

Here the problem is solved by using a high-order discretization. Having a p-
element definition automatically switches to a special formulation (sets 
"Cartesian Formulation = Logical True"). The command "Dead Loads = Logical 
False" is essential so that the load is defined to depend on the current 
deformation.

TO DO: The computed results are in good agreement, but the nonlinear iteration
denies to converge for the high loads corresponding to steps > 15 and thus 
needs some remedy.

