Elmer Release Notes for version 26.1
====================================

Previous release: **9.0**  
Period covered: **Nov 11, 2020 - Jan 19, 2026**  

These release notes provide information only on the most essential changes. Over the period there have been ~3500 commits (merge commits excluded). 
You can get a complete listing of commit messages, for example, with the command git log --since="2020-11-11"  > log.txt

It is unfortunate that the new release has been lagging. 
The detail level of release notes had to be dropped to get this out. We hope to get back on the phase.
If you want to enhance the release notes feel free to update this file and make a pull request. 

Apart from the core Elmer team at CSC (Juhani K., Mika M., Juha R., Peter R., Thomas Z.) git log shows contributions from Markus Mützel, Saeki T., Fabien G-C, Eelis T., Rich B., Iain W., Matias Z., Rupert G., Julien B., Samuel C., Luz P., Benjamin R., Ladislav M., Monher C., Olivier G., Juris V., Joe T., Bartos Z., E. Albiter, Sergey Al., Cyrille C., Jonas T., Tuomas M., Sami I., Kevin T., Sebastian G., Jonathan V., Alihossein S., Fernando B., Lucas B., Saint W., Sami R., Kasra D., Alessandro G., Andy G., Arno M., Evangelos V., Fredrik R. to this release.    

Additionally there are many ongoing developments in several branches
that have not been merged to this release and are not therefore covered here. 

The contributions of all developers are gratefully acknowledged! The work of Markus Mützel, Saeki Takayuki and
Rich Bayless is particularly appreciated as they have been consistently contributing in the true spirit of open source. Also the humongous
merge coming from univ. of St. Andrews mainly by Iain Wheel and Joe Todd related to remeshing strategies around calving is greatly appreciated.
The features have opened the door for more extensive use of adaptivity and remeshing in the code. 


### New Versioning Scheme

From this version onward we migrate to Calendar Versioning such that
- The first number (major) is the year of the 21st century, e.g. 26
- The second number (minor) is an ordinal number of releases in that year
- The third number (micro) is a growing number which for releases is always 0 and may be omitted.  


I. New Solver Modules
------------------

### FilmFlowSolver.F90
- Module for the solution of reduced dimensional Navier-Stokes equation.
- It is assumed that the average velocity of the film defines the full velocity profile and hence the depth direction may be eliminated from the flow solution.
- The intended use is film/channel flow but also cylindrical pipe flow is implemented.
- Includes lamina, Darcy and Manning models to define the relationship between pressure gradient and flux. 
- This fills the gap between full Navier-Stokes solver and reduced dimensional Reynolds solver.
- Intended use case is water flow under ice sheets. 
- See Models Manual for more details. 

### HydrostaticNSVec.F90
- Module for the solution of Stokes equation with a hydrostatic 1st order assumption a.k.a. Blatter-Patyn equation. 
- Utilizes multithreading and vectorization features in the same manner as IncompressibleNSVec, and inherits most keywords and physics from that solver. 
- For the equations look at (5.70) and (5.71) in Ralf Greve & Heinz Blatter: Dynamics of Ice Sheets and Glaciers, Springer, 2009. 
- Intended use case is 3D simulation of large ice sheets when some features of full Stokes are ignored. 

### BatterySolver/*.F90
- Solves equations for a Lithium-Ion battery from two equations for electric potential and two equations for ion concentrations.
- These equations are strongly coupled by the Butler-Volmer equation that provides the fluxes between the phases.
- The solid phase has additionally a 1D finite element equation embedded in each finite element node of the finite element mesh.
- The equation has been used in the Master's Thesis of Timo Uinonen from Univ. of Vaasa and may be considered partly experimental.
- The robustness of the solver is not the greatest requiring quite a bit of under-relaxation.
- See Models Manual for more details. 

### CyclicConvergence.F90
- This is an auxiliary solver that can be used to study whether transient and cyclic simulation has converged. It is assumed that all the cycles are saved and therefore the values from the current cycle can be compared to the previous cycle. The initial application for this solver was the case of synchronous electrical machines. 

### EMPort.F90
- Module for computing eigen modes from a special wave equation model posed over a 2-D region, typically corresponding to an electromagnetic port
- This can be used in creating boundary conditions for 3-D models based on the vector Helmholtz equation
- Not yet documented in Elmer Models Manual

### NodeToEdgeSolver.F90
- A solver that may be used to map results from nodes to edges.
- Intended use case is 3D simulation of electrical machines when the initial guess has been solved using nodal 2D approximation. 


### StatElecSolveVec.F90
- Vectorized version of StatEleSolve.

### SunAngle.F90
- Simple solver that for each node computes the maximum elevation at which the sun is still seen from a given direction.
- Intended use is related to ice sheets

### CahnHilliard.F90
- New equation for Cahn-Hilliard interface equations (with minimal testing).

### TopoOpt.F90
- Runner for topology optimization.
- Includes many filter types of which the PDE filter works well also in parallel and enables large parallel runs. 
- Has been tested in conjunction with elasticity and heat equation but the implementation should be generic.

### VectorHelmholtzNodal.F90
- Creates a discretization of the vector Helmholtz equation by using nodal (H1-conforming) finite elements.
- It is intended for creating preconditioners to solve the same equation with curl-conforming (edge) finite elements.

### WVectorFix.F90

- This contain subroutines which help in giving current sources for closed coils.
  To avoid the problem of creating a discontinuous potential, two different potential fields are solved (say
  a left and a right field). The two fields have different boundary conditions that are automatically 
  set on the bulk nodes so that a current is induced and the union of the solutions is reasonable.
  The main assumption is that the coil axis is aligned with the z-axis. 


II. Enhancements of existing solvers
--------------------------------

### CircuitsAndDynamics.F90

- A new command "Coil Use W Vector = True" to give a wire density vector
to model a stranded coil in transient/harmonic analysis (see the test case circuits_transient_stranded_wvector).
The W vector variable can be named with the command "W Vector Variable Name".
- A new command "Export Circuit Variables = Logical True" to save circuit variables in a way that
enables functional dependencies on them
- A component type "resistor" added for transient circuits
- The circuit driven solution based on the component voltage and current together with
  "Layer Electric Conductivity" and  "Layer Relative Permeability" (given in a boundary condition section).
  See for example the test case .../tests/mgdyn_harmonic_wire_impedanceBC_circuit/IBC_circuit.sif

### CoordinateTransform.F90

New keywords for creating a local coordinate system without the direction solvers:
- In a body section "Local Coordinate System Beta Reference and Gamma = Logical True"
  so that CoordinateTransform computes alpha = beta x gamma
  where beta is defined by "Beta Reference" and
  gamma is defined by "coilcurrent e" computed by CoilSolver.
- In a component section "Coil Normal(3) = Real 0.0 0.0 1.0"
  for informing CoilSolver how the coil is oriented
- In a component section "Desired Current Density = Real 1.0"
  so that CoilSolver will set a certain current density over the coil
- In a component section "Coil Use W Vector = Logical True"
  so that a coil uses a given W vector whose name can be specified as
  "W Vector Variable Name = String CoilCurrent e".
  The W vector set here is used for driving the component with circuits.
- In a boundary section "Coil Start = Logical True" as
  CoilSolver needs to know where the coil starts (if not closed)
- In a boundary section "Coil End = Logical True" as
  CoilSolver needs to know where the coil ends (if not closed)


### EMWaveSolver.F90
- The computation of eigenmodes enabled 

### FindOptimum.F90
- Nonlinear optimization methods have been added.
- Elmer can control itself when performing optimization, and does not need external optimization package. 


### GmshOutputReader.F90
- This reader can read ascii output written in Gmsh format.
- The original use case was to read data on boundaries that hierarchically couple to other BCs.

### HeatSolveVec.F90
- Discontinuities enabled between bodies.
- Enable modeling of radiation heat transfer using the concept of "Radiosity". 
- Enable computation of diffuse gray radiative heat transfer in parallel.
- Except for the 1st feature these are available also in the legacy HeatSolve.F90 

### HelmholtzSolve.F90
- Preconditioning of shifted Laplacian type is now possible

### MagnetoDynamics
  - A-V solution with a user-defined reluctivity function
  - A model of Darwin type enabled
  - A new command "Calculate Homogenization Loss" for stranded coils modelled in the frequency domain
  - An option to convert a H-B curve to effective values for harmonic analysis as
  
    H-B Curve = Variable time ! ... or anything else always available (this is not used in the evaluation)  
       Real Harmonic ! cubic monotone ...  
         include HB  
       End

### MagnetoDynamics2D.F90
  - The implementation of the London equations added; see the test cases fem/tests/circuits2D_*_london

### ParticleAdvector.F90
  - Improved robustness and path integration 

### ResultOutputSolve
  - High-order visualization for the p-version of FEM
  - Gmsh alternative to the default file I/O based coupling strategy
  - First version of saving surface mesh in STL format. 


### SaveData
  - Enable saving of data on segments of lines also with p- and H(curl) elements.

### SaveGridData.F90
  - Support for NETCDF output files

### ShellSolver.F90
  - Full geometric nonlinearity with low-order shell elements
  - Higher-order shell elements backgrounded by a standard surface mesh

### Smitc.F90
  - A 4-node plate element based on the full linked interpolation is also available

### VectorHelmholtz.F90
  - Block preconditioning based on the hierarchic grouping of DOFs (for the quadratic approximation)
  - Block preconditioning with constraints enables
  - "Mass-Proportional Damping = Logical True" can be used to create preconditioners based on a mass-proportional perturbation 
  - Alternate formulations based on using a scalar potential have been implemented
  - 2D simulations enabled
  - The second-order Nedelec bases of the 2nd kind for tetrahedra and triangles enabled
  - Ability to create ports automatically.
  


III. ElmerSolver library functionality
---------------------------------

### Linear solvers
  - Additional mumps solvers beyond dmumps, i.e., zmumps, cmumps, smumps
  - AMGX interface can be used by giving "Linear System Solver = amgx"
  - An interface to Rocalution library may be used by giving "Linear System Solver = rocalution"
  - Possibility to use Hypre library also with one MPI process
  - p-multigrid for H(curl) discretizations of order 2. Currently this works best for
    the second-order Nedelec bases of the 2nd kind (available for tetrahedra and triangles)
  - Another solver can be used for the preconditioning by using the keyword "Prec Solvers".
    The preconditioning solver must be able to read as input the residual variable (whose name is here specified with the keyword
    "Preconditioning Residual") and produce the correction variable (whose name is here specified as the value of the
    keyword "Preconditioning Update").
  - "Linear System Preconditioning = Cholesky" as a synonym for
    "Linear System Symmetric=True" and "Linear System Symmetric ILU=True"
  - The normwise relative backward error err = ||r||/(||A|| ||x|| + ||b||) may be used as
    a convergence criterion also for complex linear systems
  - New smoothers in conjunction with multigrid and multilevel methods. 

### Solvers for eigenproblems
  - "Linear System Direct Method = cholmod" allowed in the solution of eigenproblems

### Time integration
  - "Generalized-alpha" time integration method for systems with the time derivatives of order 2
     may be activated by "Timestepping Method = Generalized-alpha".
    The spectral radius to define "optimal" parameters may be given
    with "Generalized-alpha rinf = Real ..." (0-1, 0.75 as default).

### Finite elements
   - p-element basis functions based on complete polynomials can now be used by giving
     "Serendipity P Elements = False" in the simulation section. The default value "True" leads to using
     serendipity bases.
  - A 14-node pyramid element added
  - The second-order Nedelec bases of the 2nd kind for tetrahedra and triangles to enable curl-conforming approximation with
    12-DOF triangle and 30-DOF tetrahedron. Note that similar approximations for other element shapes have not yet been
    implemented.
  - Basis functions of a 48-DOF div-conforming vector element for bricks
  - Div-conforming basis functions for higher-order Raviart-Thomas RT_1 (8-DOF) element over triangles 
  - An element definition "n:N ...", N>1 allowed, although some features are not fully supported.
    BCs may be given by using a part {n} in the BC definition.
    For example, if we have
    
       Variable = E  
       Variable DOFs = 1  
       Element = "n:1 ..."
    
    then a BC definition
    
       E {n} = Real ...
    
    should also work. On the other hand, if
    
       Variable = E  
       Variable DOFs = 1  
       Element = "n:2 ..."
    
    then we may give
    
       E {n} 1 = Real ...  
       E {n} 2 = Real ...
    
    Finally, if we the variable is defined to have components and the component
    names are created for example from
    
       Variable = E[E Re:1 E Im:1]  
       Element = "n:2 ..."
    
    we may define
    
      E re {n} 1 = Real ...  
      E re {n} 2 = Real ...  
      E im {n} 1 = Real ...  
      E im {n} 2 = Real ...  

    Note that this is still simplistic implementation and lacks many options which work
    in the case of standard BCs for nodal (Lagrange) interpolation approximations.
  - Matrix representations of Nedelec's interpolation operators added


### Constraints
  - Automated naming of Lagrange multipliers so that each solver has its own name. The command "Lagrange Multiplier Name"
    can also be used to give the name.
  - Functional changes related to mortar keywords: the keywords "Galerkin Projector" or "Level Projector Strong"
    should be used to select whether the mortar constraints are generated weakly or not.
    If this information is not provided, the constraints are now generated in a strong manner corresponding to
    the choices "Galerkin Projector = False" and/or "Level Projector Strong = True".
    For a more fine-tuned selection the keyword "Level Projector Nodes Strong" is
    also available. On the other hand the keyword "Level Projector Generic = True" can be used to
    select the newest subroutines for generating mortar constraints. This has effects
    on the generation of both weak and strong constraints.
  - MPRGP iterative method for problems with non-equality constraints

### Optimization
  - New external optimization routines: a square finding routine from minpack, Bobyqa and Newuoa
    from Powell's optimization routines. The code was modified so that the optimization
    methods can basically call Elmer to provide the cost function - not the other way round. This requires
    that the basic solution cycle is provided as a pointer to a function. See also the compilation
    flag "WITH_EXTOPTIM".

### Miscellaneous
  - Option to use MMG* library internally for adaptive mesh generation
  - Internal extrusion into separate MPI tasks.
  - Internal splitting of prismatic meshes into tetrahedrons. 
  - Tentative machinery for CutFEM without elimination implemented
  - Possibility to create BC's on-the-fly based on geometric detection.
  - Feature to "calculate mesh pieces" to ensure conformity of meshes. 


IV. ElmerGrid
-------------
- Add reading Gmsh input format 4.1 in binary.  
- Add calculate mesh pieces function, which checks for conformity of each mesh.
- Report size of bounding box, and give hint to user if mesh may be in millimeters.
- Numerous small fixes.

V. ElmerGUI
-----------
- Users can add their own material libraries by putting xml files in edf-extra directory.
- Tango-based icons for better visibility in high resolution monitors
- Add reading Gmsh input format 4.1 in binary.
- Add reading STL input format in binary.
- Numerous small improvement and fixes
- Library versions updated


VI. Configuration & Compilation
--------------------------------
- A new flag "WITH_MMG" to compile with MMG
- A new flag "WITH_EXTOPTIM" to include external optimization routines
- A new flag "WITH_AMGX" to use the AMGX interface
- A new flag "WITH_ROCALUTION" to compile with Rocalution library
- A CHOLMOD option by having -DWITH_CHOLMOD:BOOL=TRUE and by providing
  CHOLMOD_LIBRARIES and CHOLMOD_INCLUDE_DIR
- An option to activate an external UMFPACK by using EXTERNAL_UMFPACK, UMFPACK_LIBRARIES and UMFPACK_INCLUDE_DIR
- USE_SYSTEM_ZOLTAN to pick external or compiled library
- A new flag "HAVE_QP" for indicating quadratic-precision real numbers


VII. Elmer/Ice
--------------
- New features in Elmer/Ice are documented elsewhere


VIII. Obsolete code
------------------
- ElmerPost moved to separate repository "ElmerPost"
- Obsolete Trilinos interface removed
