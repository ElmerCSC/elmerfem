Elmer Release Notes for version 8.5
===================================

Previous release: **8.4**  
Period covered: **18 Dec 2018 - 30 Aug 2020**  
Number of commits: **~1140** (excluding merges)  

These release notes provide information on the most essential changes.
You can get a complete listing of commit messages, for example, with:  
git log --since="2018-12-18"  > log.txt

Apart from the core Elmer team at CSC (Juhani K., Mika M., Juha R., Peter R., Thomas Z.)
git log shows contributions from from Daniel B., Denis C., Eef v. D., Eelis T., Fabien G-C,
Foad S. F., Fredrik R., Olivier G., Joe T., Luz P., Mondher C., Rupert G., Sami I.,
Sami R., Samuel C., and Saeki T. have contributed to making this release. 

Additionally there are many ongoing developments in several branches
that have not been merged to this release, and are therefore not covered here. 
Also sometimes the code has been passed on by the original author by means other than the
git, and in such cases the names may have been accidentally omitted.

The contribution of all developers is gratefully acknowledged! 


New Solver Modules
------------------

### IncompressibleNSVec
- Incompressible Navier-Stokes solver utilizing vectorized and threaded assembly
- Includes built-in support for block preconditioning (Schur complement approximation included)
- Includes non-Newtonian material laws
- Intended for Elmer/Ice community but also other may find it useful. 


### BeamSolver3D
- Solver for the Timoshenko equations of elastic beams embedded in 3-D space (see Elmer Models Manual for documentation)

### GmshReader
- Reads the mesh and results from simple Gmsh file format (that can be written by ElmerSolver as well)
- Solver includes interpolation of the fields to the current mesh
- May be used for hierarchical simulations where results are inherited from previous simulations

### ModelMixedPoisson
- A general-purpose mixed FEM solver for the Poisson equation (see Elmer Models Manual for documentation)
- employs a div-conforming (face) finite element approximation

### SpringAssembly
- A generic utility to add node-wise springs and masses to structural models (see Elmer Models Manual for documentation)


Enhanced Solver Modules
-----------------------

### ElasticSolve
- adding a new UMAT material model is simplified: compilation with an elmerf90 command is sufficient
- the state variables of UMAT material model can be written to a result file and visualized
- UMAT implementation updated to support axial symmetry 

### EMWaveSolver
- the solver updated to support the basis functions of second order and simulation in 2D 
- the solver is now documented in Elmer Models Manual

### MagnetoDynamics2D
- a velocity field can be given to add a Lorentz term to the equations

### ResultOutputSolver
- Vtu format:
    - Enable saving of pieces, i.e. bodies and boundaries 
    - Improved saving of elemental, DG and IP fields
- Gmsh format:
    - Improved use of masking features in output

### ShellSolver
- Eigenanalysis with the shell solver enabled
- Spring, resultant force and couple BCs added
- Combined analysis of 2-D shells and 1-D beams enabled
- Fully coupled analysis of 2-D shells and 3-D solids enabled (still subject to some geometric constraints on the mesh) 

### StructuredMeshMapper
- Enable arbitrary number of layers, before limited to three. 

### HeatSolver
- A new tentative vectorized version: HeatSolverVec
- Enable symmetric 3D cases for view factor computation. 
- Make Gebhart factors linear system symmetric, if possible "ViewFactor Symmetry".

### StressSolver
- Added Maxwell visco-elastic model to linear elasticity solver
- possible also to be run incompressible (introducing pressure variable)
- optional pre-stress advection term for layered Earth-deformation model

### WhitneyAVHarmonicSolver
- A surface impedance condition for the time-harmonic AV model

ElmerSolver library functionality
---------------------------------

### Treatment of block systems
- The block matrix approach for solving complicated problems has been enhanced.
  Currently the block approach can be used in several ways during some stage of the solution.
     1. Split up monolithic equations into subproblems that are easier to solve (e.g. IncompressibleNS)
     2. Combine linear multiphysical (coupled) problems into a block matrix (e.g. FSI problems)
- For problems belonging to class 1) we may perform recreation of a monolithic matrix. This will
  allow better use of standard linear algebra to utilize direct solvers, or change the system to
  be harmonic or eigenvalue problem.
- For the documentation of utilizing block-matrix construct in connection with
  the fully coupled simulation of multiphysical problems see the new chapter
  "Block-matrix construct to build tightly coupled solvers" in ElmerSolver Manual. 

### More economical integration rules
- A collection of economical Gauss quadrature rules for prismatic elements are introduced to replace
    tensor product rules for quadrilateral p-elements when 1 < p <= 8. The tensor
    product rule with n = (p+1)**2 points is now replaced by more economical ones.

### Dirichlet BCs for div-conforming vector finite elements (face elements)
-  A sif command of the form Q {f} j = Real ... can be used to specify vector-valued data whose 
     normal component is then used to integrate the values of DOFs for vector-valued interpolation of the data.
     Here Q is an Elmer variable which is approximated with face finite elements.

### Conforming BCs by elimination
- System can identify conforming boundaries such that dofs related to nodes or edges on opposing sides may be
  assembled into one degree of freedom.
- This decreases the size of the linear system and is numerically favourable.
- Antiperiodicity may be included. For vector-valued problems all components must be treated alike. 
- Conforming BCs for edge dofs may consider direction of edge.
- See test cases with "Apply Conforming BCs" and "Conforming BC" defined.
    
### Improved internal partitioning with Zoltan
- Enable internal partitioning with Zoltan to honor connected boundaries.

### Enable primary solver to call other solvers
- For documentation see the section "Solver execution by a master solver" in ElmerSolver Manual.
- Enables calling before and after solving the primary problem.
- Also possible to call before and after each nonlinear iteration.


### Anderson Acceleration for linear systems
- Implemented a version of Anderson Acceleration where previous solutions and
  residuals are used to accelerate the nonlinear convergence.
- May increase linear convergence to quadratic, quadratic convergence (Newton's method) is not improved.

### Run Control
- Enable external loop control over the simulation.
- May be used in optimization and parametric scanning etc.
- Applicable also to transient system as the variable "time" is not used for the control level. 



ElmerGrid
---------
- Fixes for UNV, mptxt and Gmsh file format import.
- Tentative reader for FVCOM format
- Add possibility to define seed for Metis partitioning (-metisseed).
- Maintain entity names in extrusion



ElmerGUI
--------
- Include many improvements by Saeki

Configuration & Compilation
---------------------------
    
- cmake improvements:


Elmer/Ice
---------
- New features in Elmer/Ice are documented in elmerfem/elmerice/ReleaseNotes/release_elmerice_8.5.md



