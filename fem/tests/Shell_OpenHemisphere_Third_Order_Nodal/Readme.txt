A benchmark problem where an open hemispherical shell is subject to point
loads. The reference results can be found in Sze KY, Liu XH, Lo SH. Popular 
benchmark problems for geometric nonlinear analysis of shells. Finite 
Elements in Analysis and Design 2004, 40(11):1551-1569.

Here the problem is solved by using a third-order nodal mesh of the physical 
shell. In this case all geometric information is derived from the mesh and no
user-supplied information about the shell director is needed. At the moment 
no numerical tricks are applied to handle numerical over-stiffness (locking),
but the basic third-order approximation may give reasonable results if the 
shell is not very thin. This approach needs a command "Mesh Reparametrization =
Logical True".

For approximation with p-elements see the test 
Shell_OpenHemisphere_High_Order_Blending
