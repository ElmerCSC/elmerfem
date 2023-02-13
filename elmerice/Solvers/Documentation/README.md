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

# Providing documentation for ElmerIce Solvers

Github already renders .md files but these do not support maths very well.  
You can generate a .pdf that will display the maths better (if you've written it in the .md file using LaTEX syntax) using e.g. pandoc.
> pandoc -f gfm -t pdf -s SomeDocumentationFile.md > SomeDocumentationFile.pdf
will generate a fully rendered pdf version of any of the files.

If you've written a new solver, please use the [Template.md](./Template.md) file to set up your documentation before uploading it to this directory.

To generate the doucumentation for the adjoint method run:
> pandoc -d MakeDoc_Adjoint.yml


## Adjoint methods

Following a re-organisation of the adjoint methods in Elmer/Ice in April 2020,  

- generic solvers are found under [elmerice/Solvers/Adjoint](../Adjoint)
- model dependent solvers are found under elmerice/Solvers/Adjoint*[Name_Of_Model]*,
where *[Name_Of_Model]* can be Stokes, SSA, Thickness (to come...).

A table showing the old solvers and their replacement can be found 
[here](https://cloud.univ-grenoble-alpes.fr/index.php/s/AHCwsgKgjWimqdG).
A fatal message is now called in *old* solvers to advertise for the replacements and these solvers will be removed in future releases.

Here is the list of supported and documented solvers:
 
Generic Solvers:  

- [Optimize_m1qn3](Optimize_m1qn3.md)   
- [Adjoint_GradientValidation](Adjoint_GradientValidation.md)
- [Adjoint_LinearSolver](Adjoint_LinearSolver.md)
- [Adjoint_CostDiscSolver](Adjoint_CostDiscSolver.md)
- [Adjoint_CostContSolver](Adjoint_CostContSolver.md)
- [Adjoint_CostRegSolver](Adjoint_CostRegSolver.md)

Stokes Solvers: 

- [AdjointStokes_GradientBetaSolver](AdjointStokes_GradientBetaSolver.md)

SSA Solvers:

- [AdjointSSA_SSASolver](AdjointSSA_SSASolver.md)
- [AdjointSSA_GradientSolver](AdjointSSA_GradientSolver.md)
- [AdjointSSA_CostFluxDivSolver](AdjointSSA_CostFluxDivSolver.md)
- [AdjointSSA_CostTaubSolver](AdjointSSA_CostTaubSolver.md)

Mass Conservation (thickness solver):

- [AdjointThickness_ThicknessSolver](AdjointThickness_ThicknessSolver.md)
- [AdjointThickness_GradientSolver](AdjointThickness_GradientSolver.md)


Generic user functions:

- [Utility](Utility.md)   

## Coupled hydrology-plumes-calving

This is pretty much the output of Samuel Cook's PhD thesis on 3D coupled modelling of a tidewater glacier, and involves coupling the GlaDS hydrology solvers with Joe Todd's 3D calving solvers, and a new 1D ODE solver for glacial meltwater plumes based on the work of Donald Slater. Several other new solvers are also required to manage the interaction between all these moving parts. If you're interested in using this set-up, a full description, including all necessary solvers, SIF inclusions, and mesh fiddliness, is provided in [CoupledIceHydrologyCalvingPlumesDocumentation](CoupledIceHydrologyCalvingPlumesDocumentation.md)(the individual new solvers are also all documented in their own .md files, which you may find it useful to look at). Note: all the necessary modifications to existing Elmer/Ice solvers are in the distributed versions, so you shouldn't have to do anything not listed in the doc, as long as I've not forgotten to tell you about something important. Any questions, email me at samuel.cook .at. univ-grenoble-alpes.fr
