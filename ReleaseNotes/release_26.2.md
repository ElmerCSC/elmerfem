Elmer Release Notes for version 26.2
====================================

Previous release: **26.1**  
Period covered: **Jan 19, 2026 - April 24, 2026**  

These release notes provide information only on the essential changes.
Over the period there has been ~117 commits (merge commits excluded). 
Some small fixes, less relevant changes and plain code refactoring have been omitted.
With this release we are back in pace where it is easier to summarize the recent changes. 
You can get a complete listing of commit messages, for example, with the command git log --since="2026-19-1"  > log.txt

Apart from the core Elmer team at CSC (Juhani K., Mika M., Juha R., Peter R., Thomas Z.) git log shows contributions from Markus Mützel,
Fabien G-C, Rupert G.,C. Bosbeux, K. Teenakul, B. Rodenberg, D. Bates and T. Mylläri (including changes to build system).  
Additionally there are many ongoing developments in several branches
that have not been merged to this release and are not therefore covered here. 

The contributions of all developers are gratefully acknowledged!


I. New solver modules
---------------------

### ADIOS2OutputSolver
- Saves elmer nodal scalar and vector fields in ADIOS2/BP5 format and outputs Fides json file for Paraview.
  This greatly reduces the number of files in parallel runs. 
 

I. Enhancements of existing solvers
-----------------------------------
 
### EMPort.F90
- Create a separate routine for just solving the port potential that inherits the material parameters from the parent material.
- Make element definitions to follow primary solver
- Have port conditions for vector Helmholtz equation only in one place.
- Enable eigenmodes to be sorted after solution with "Eigen System Sorting = String". Same values available as for "Eigen System Select".
- In EMPortSolver enable solution of several ports such that one port is solved at a time. Currently probably limitations in parallel.
- Add "Port Type" value "eigenmode" that automatically picks the eigenmodes computed by EMPortSolver.
- Add keyword "Port Passive" that turns off the r.h.s. of BC.
- When doing "Eigen System Scale to Unity" divide with the max abs complex number such that also the eigenmode is rotated at the same time.
- In EMPortSolver enable automated shifting based on analytic limit of the propagation constant. Also save the propagation constant of each boundary for later use.

### MarchingODE.F90
- subdivision of timestep in marching ode solver for improved accuracy
- internal variables to be used in marching.
- vector valued field in marching
- Enable MarchingODESolver to have multiple components.

### CalcFields.F90
- Optionally allow Lorenz force to be computed as a nodal force instead of being a distributed one.



III. ElmerSolver library functionality
---------------------------------------

### Finite elements

- Test 2nd-order prisms and p-multigrid preconditioning for H(curl)
-Alternate H(curl) basis functions of degree 2 for prisms

Some changes related to the basis functions in H(curl) (#765)    
- some restructuring to avoid code repetition
- an alternative 2nd order quad added
- some preparations to allow for the cubic approximation in the case of tetrahedra
- a 24-point quadrature for tetrahedra added

 quadratic approximation in H(curl) (#798)
    
    This relates to cases where the mesh combines several element types. Some basis
    functions associated with quad faces were not compatible over faces shared by
    elements of different types.

    Modernize calls to get the vector element basis functions

    A revision of H(curl) finite elements (#750)
    
    * The original code for constructing the lowest-order H(curl) basis functions of the second kind is replaced by the construction where some of the basis functions are obtained as the gradients of scalar fields.
    
    * Alternate versions of the second-order Nedelec bases of the first kind for simplicial elements. The motivation is that better iterative solvers can be obtained for these H(curl) approximations. A special keyword "Simplicial Mesh = Logical True" must be given to activate their use, so the appropriate combination of keywords is then as follows
    
      Second Kind Basis = False
      Quadratic Approximation = True
      Simplicial Mesh = Logical True
    
    Analogous versions for other element shapes do not yet exist.
    
    * Some corrections for setting non-homogeneous Dirichlet constraints for H(curl) approximations. A new 3-D test case related to such BCs is added. A solver employed in the testing is also moved to the modules directory.


- Fixes for elementtype 715 and added support for elementtype 718.





### Miscallenous

- Enable reading of binary mesh format written by ElmerGrid. The changes are implemented by adding a suffix 'bin' to express the binary nature. Also implements single precison saving of nodes with suffix 'sbin'. To activate these use flag '-bin' and '-sbin' in ElmerGrid. ElmerGUI not yet supported nor is the reading in ElmerGrid.
- Accurate integration over step functions using temporal splitting of elements, works only in 2D.
- Enable Hcurl Dirichlet BC's such that even if it does not share the full boundary the values can still be set at the edges.
- Function for interpolating many curves at the same time.
- Also add possibility for namespace for GaussPointsAdapt such that we can choose "bc gauss" or "bulk gauss" rules in assembly. Add test case with exact known solution.
- Display current branch in ElmerSolver
- Set the default value of global bubbles to False thereby removing the strange problems when bubbles are inherited to solvers that do not need them.
- ADIOS2Utils module with a class for simplified writing of parallel data.
- First steps to implement generic Nitsche type of boundary conditions for nonconforming interfaces. This includes tentative implementation for 1D boundaries and a related test case.



IV. Bug Fixes
--------------

- Eliminate empty space from timesteps in XML because VTK is more picky lately.
- Fix in passive parallel interfaces that appears when passive interface and partition interface overlap.
- Fix adaptive integration with temporal triangles.
- Fix segfault in mortar projectors with discontinuous Galerkin-
- Fix the restart when there is a combination of different permutation vectors and field that is not everywhere defined.
- Set locale "C" before formatting strings in MATC
- Fix initialization of halo nodes/elements when working with halo meshes
- Fix 1D weak projector in case of DG elements
 

V. ElmerGrid
-------------
- Enable binary writing of the mesh (except for the header and names)
- Display compilation date and current branch in ElmerGrid
- Support for elementtype 718 (quadratic prism with centernodes). 


VI. ElmerGUI
-----------
- Add tentative support also for elementtype 718.


VI. Configuration & Compilation
--------------------------------
- Many improvements on continious integration and testing.


VII. Elmer/Ice
--------------
- New features in Elmer/Ice are documented elsewhere



   Merge pull request #801 from ElmerCSC/elmerice    
    Elmerice


    Author: fgillet <fabien.gillet-chaulet@univ-grenoble-alpes.fr>
    ComputeGroundingLineFlux    
    Add solver to compute grounding line flux in 3D extruded simulations

       Author: RupertGladstone <rupertgladstone1972@gmail.com>
    Minor bug fix to Flotation to allow it to run on the lower surface of a 3D body.

    Cyrille:
    Add initialization of GroundedMask for all active nodes    
    Initialize GroundedMask for active DOFs to avoid default values in halo nodes, leaving wrong GroundedMask value at the partition interfaces with halo elements/nodes.


    Author: Cyrille <91067824+cmosbeux@users.noreply.github.com>
    Initialize LimitedSolution to false    
    Fix uninitialized LimitedSolution on halo elements: LimitedSolution was only visited for active elements. It should be initialised to True in passive/halos elements to avoid limiter inconsistencies and free-surface artifacts at partition boundaries.






Author: fgillet <fabien.gillet-chaulet@univ-grenoble-alpes.fr>
    Fix Dirichlet for passive Elements
    
    - PassPerm was not allocated and initialised for cases where Variable % DOFs > 1; In this case the passive mechanism is set by "VarName Passive = Logical True" but Dirichlet conditions are set for each component individually.
    - Do not set Dirichlet conditions from halo elements. Correct dettection of the passive/active boundary requires halo elements; but might happen that it is not detected for halos at the border, resulting in wrong Dirichlet conditions.


    Author: Juha Ruokolainen <jpr@keisarikotka.lnx.csc.fi>
    - disable "ForceToStress_parallel" elmerice test if no parallel direct solvers
    - fix uninitialized access in elmerice "ComputeNormal()" solver
    - disable (very old) "EliminateDirichlet" in "contact" test, doesn't interact
    well with "ComputeLoads"


IX. Obsolete code
------------------
- Removed some mortar projectors assuming the "stride" strategy
