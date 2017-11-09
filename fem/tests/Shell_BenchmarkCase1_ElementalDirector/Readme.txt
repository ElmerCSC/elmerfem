Use mesh.elements.data to define the shell director data.

! Solve nonlinear shell equations of Reissner-Naghdi type in the case of 
! a cylindrical benchmark problem with membrane-dominated asymptotic behaviour in 
! the linear regime. The linearized problem has been described in 
! Pitkäranta et al. Shell deformation states and the finite element method: 
! a benchmark study of cylindrical shells. Computer Methods in Applied Mechanics 
! and Engineering 1995. 128:81-121. As opposed to the reference, here the shell
! problem is treated as nonlinear but the shell is subject to a small load 
! so that nonlinear effects are not yet significant.

Results below: The error of strain energy when the shell problem is linear
(Large Deflection = False) and when the mesh is refined as


Element Divisions 1 = 1 16 1  
Element Divisions 2 = 1 8 8 1

I. QUADS

Strain reduction method = 1 (solver default)
================================
8.5183751273725697E-002
4.7861369052476552E-002
2.6866602832960274E-002
1.5742186303373210E-002
9.7455027946871976E-003

Strain reduction method = 0 (no reduction)
======================================
0.51900858861935217
0.30075442699998356
0.16281741500325506
8.4638690386812335E-002
4.3193057630034115E-002
