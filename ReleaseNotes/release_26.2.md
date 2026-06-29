Elmer Release Notes for version 26.2
====================================

Previous release: **26.1**  
Period covered: **Jan 19, 2026 - April 28, 2026**  

These release notes provide information only on the essential changes.
Over the period there have been ~128 commits (merge commits excluded). 
Some small fixes, less relevant changes and plain code refactoring have been omitted.
With this release we are back in pace where it is easier to summarize the recent changes. 
You can get a complete listing of commit messages, for example, with the command git log --since="2026-19-1"  > log.txt

Apart from the core Elmer team at CSC (Juhani K., Mika M., Juha R., Peter R., Thomas Z.) git log shows contributions from M. Mützel,
Fabien G-C, Rupert G.,Cyrille M., Kavin T., Benjamin R., Daniel B. and Tuomas M. (including changes to build system).  
You can check the authors related to this release, for example, with (note that squashes may destroy author contributions), 
git shortlog --since="2026-19-1" -nes

Beyond this release there are many ongoing developments in several branches (e.g. on constraint solutions, radiation heat transfer)
that have not been merged to this release and are not therefore covered here. 

The contributions of all developers are gratefully acknowledged!


I. New solver modules
---------------------

### ADIOS2OutputSolver
- Saves Elmer nodal scalar and vector fields in ADIOS2/BP5 format and outputs Fides json file for Paraview.
  This greatly reduces the number of files in parallel and time-stepping runs.
 

II. Enhancements of existing solvers
-----------------------------------
 
### EMPort.F90
- Create a separate routine for just solving the port potential that inherits the material parameters from the parent material.
- Make element definitions to follow the primary solver
- Have port conditions for vector Helmholtz equation only in one place.
- Enable eigenmodes to be sorted after solution with "Eigen System Sorting = String". Same values available as for "Eigen System Select".
- In EMPortSolver enable solution of several ports such that one port is solved at a time. Currently probably limitations in parallel.
- Add "Port Type" value "eigenmode" that automatically picks the eigenmodes computed by EMPortSolver.
- Add keyword "Port Passive" that turns off the r.h.s. of BC.
- When doing "Eigen System Scale to Unity" divide with the max abs complex number such that also the eigenmode is rotated at the same time.
- In EMPortSolver enable automated shifting based on analytic limit of the propagation constant. Also save the propagation constant of each boundary for later use.

### MarchingODE.F90
- subdivision of timestep in marching ODE solver for improved accuracy
- internal variables to be used in marching.
- Enable MarchingODESolver to have multiple components.

### CalcFields.F90
- Optionally allow Lorenz force to be computed as a nodal force instead of being a distributed one.



III. ElmerSolver library functionality
---------------------------------------

### Finite elements

- Alternative versions of the second-order Nedelec bases of the first kind for triangles, quads, tetrahedra and prisms. The motivation is that better iterative solvers can be obtained for these H(curl) approximations. A special keyword "Gradient Basis Functions = True" must be given to activate their use. Analogous versions for other element shapes do not yet exist.
- The original code for constructing the lowest-order H(curl) basis functions of the second kind is also replaced by the alternative construction.
- The third-order H(curl) basis for triangles (the Nedelec first family).
- A 24-point quadrature for tetrahedra added.
- Fixes for the 15-node prism (element type 715) and support added for the 18-node prism (the element type 718).


### Miscellaneous

- Enable reading of binary mesh format written by ElmerGrid. The changes are implemented by adding a suffix 'bin' to express the binary nature. Also implemented is single precision saving of nodes with suffix 'sbin'. To activate these use flags '-bin' and '-sbin' in ElmerGrid. ElmerGUI not yet supported nor is the reading in ElmerGrid.
- Accurate integration over step functions using temporal splitting of elements, works only in 2D.
- Enable Dirichlet BCs for H(curl) such that even if it does not share the full boundary the values can still be set at the edges.
- Function for interpolating many curves at the same time.
- Also add possibility for namespace for GaussPointsAdapt such that we can choose "bc gauss" or "bulk gauss" rules in assembly. Add a test case with exact known solution.
- Display current branch in ElmerSolver (requires compilation from git repository). 
- Set the default value of global bubbles to False thereby removing the strange problems when bubbles are inherited to solvers that do not need them.
- ADIOS2Utils module with a class for simplified writing of parallel data.
- First steps to implement generic Nitsche type of boundary conditions for nonconforming interfaces. This includes tentative implementation for 1D boundaries and a related test case.



IV. Bug Fixes
--------------

- Eliminate empty space from timesteps in XML because VTK is more picky lately.
- Fix in passive parallel interfaces that appears when passive interface and partition interface overlap.
- Fix adaptive integration with temporal triangles.
- Fix segfault in mortar projectors in combination of discontinuous Galerkin
- Fix the restart when there is a combination of different permutation vectors and a field which is not everywhere defined.
- Set locale "C" before formatting strings in MATC
- Fix initialization of halo nodes/elements when working with halo meshes
- A fix for quadratic approximation in H(curl) relating to cases where the mesh combines several element types. Some basis functions associated with quad faces were not compatible over faces shared by elements of different types.
 

V. ElmerGrid
-------------
- Enable binary writing of the mesh (except for the header and names)
- Display compilation date and current branch in ElmerGrid
- Support for the element type 718 (quadratic prism with center nodes). 


VI. ElmerGUI
-----------
- Add tentative support also for element type 718.


VII. Configuration & Compilation
--------------------------------
- Many improvements on continuous integration and testing.


VIII. Elmer/Ice
--------------
- Add solver to compute grounding line flux in 3D extruded simulations
- Minor bug fix to Flotation to allow it to run on the lower surface of a 3D body
- Add initialization of GroundedMask for all active nodes    
- Initialize GroundedMask for active DOFs to avoid default values in halo nodes, leaving wrong GroundedMask value at the partition interfaces with halo elements/nodes
- Initialize LimitedSolution to false    
- Fix uninitialized LimitedSolution on halo elements: LimitedSolution was only visited for active elements. It should be initialised to True in passive/halos elements to avoid limiter inconsistencies and free-surface artifacts at partition boundaries
- Fix Dirichlet for passive Elements
    - PassPerm was not allocated and initialised for cases where Variable % DOFs > 1; In this case the passive mechanism is set by "VarName Passive = Logical True" but Dirichlet conditions are set for each component individually.
    - Do not set Dirichlet conditions from halo elements. Correct detection of the passive/active boundary requires halo elements; but might happen that it is not detected for halos at the border, resulting in wrong Dirichlet conditions.
- disable "ForceToStress_parallel" elmerice test if no parallel direct solvers
- fix uninitialized access in elmerice "ComputeNormal()" solver
- disable (very old) "EliminateDirichlet" in "contact" test, doesn't interact
    well with "ComputeLoads"


IX. Obsolete code
------------------
- Removed some mortar projectors assuming the "stride" strategy
