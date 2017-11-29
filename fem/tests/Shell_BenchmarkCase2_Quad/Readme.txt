Solve nonlinear shell equations of Reissner-Naghdi type in the case of 
a cylindrical benchmark problem with bending-dominated asymptotic behaviour in 
the linear regime. The linearized problem has been described in 
Pitkäranta et al. Shell deformation states and the finite element method: 
a benchmark study of cylindrical shells. Computer Methods in Applied Mechanics 
and Engineering 1995. 128:81-121. As opposed to the reference, here the shell
problem is treated as nonlinear but the shell is subject to a small load 
so that nonlinear effects are not yet significant.

Results below: The error of strain energy when the shell problem is linear
(Large Deflection = False) and when the mesh is refined as

Element Divisions 1 = 1 16 1  
Element Divisions 2 = 1 8 8 1

I. QUADS

The default (strain reduction) method
================================
0.28534440085557611
0.14033086583199897
7.0275522871128052E-002
3.5170214283552041E-002
1.7592769166581007E-002


Strain reduction method = 0 (no reduction)
======================================
0.99896892914864110
0.99574227718518038
0.98313882259801899
0.93712411256923966
0.80190890167786044
