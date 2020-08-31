# Providing documentation for ElmerIce Solvers

Github already render .md files but do not support math very well.  
You can generate .pdf using e.g. pandoc.

To generate the doucumentation of the adjoint method run:
> pandoc -d MakeDoc_Adjoint.yml

## Adjoint methods

 Following a re-organisation of the adjoint methods in Elmer/Ice in April 2020,  

- generic solvers are found under [elmerice/Solvers/Adjoint](../Adjoint)
- model dependent solvers are found under elmerice/Solvers/Adjoint*[Name_Of_Model]*,
where *[Name_Of_Model]* can be Stokes, SSA, Thickness (to come...).

A table showing the old solvers and their replacement can be found 
[here](https://cloud.univ-grenoble-alpes.fr/index.php/s/AHCwsgKgjWimqdG).
A warning message will be  dispalyed by *old* solvers to advertise for the replacements and theses solvers will be removed in future releases.

Here is the list of supported and documented solvers:
 
Generic Solvers:  

- [Optimize_m1qn3](Optimize_m1qn3.md)   
- [Adjoint_GradientValidation](Adjoint_GradientValidation.md)
- [Adjoint_LinearSolver](Adjoint_LinearSolver.md)
- [Adjoint_CostDiscSolver](Adjoint_CostDiscSolver.md)
- [Adjoint_CostContSolver](Adjoint_CostContSolver.md)
- [Adjoint_CostRegSolver](Adjoint_CostRegSolver.md)

SSA Solvers:

- [AdjointSSA_SSASolver](AdjointSSA_SSASolver.md)
- [AdjointSSA_GradientSolver](AdjointSSA_GradientSolver.md)
- [AdjointSSA_CostFluxDivSolver](AdjointSSA_CostFluxDivSolver.md)
- [AdjointSSA_CostTaubSolver](AdjointSSA_CostTaubSolver.md)

Generic user functions:

- [Utility](Utility.md)   

