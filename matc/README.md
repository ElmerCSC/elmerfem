#MATC libary 

MATC is a library for the numerical evaluation of mathematical expressions.
In the Elmer package it can be used in defining expressions in the command
file of ElmerSolver. The expressions may be evaluated either while the
command file is read in, or at the time of execution.

MATC is ideally used for quick testing or things that are not done too
often. For performance critical applications rather use UDF's
(user defined function) or even Lua. 