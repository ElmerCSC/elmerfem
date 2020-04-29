# Providing documentation for ElmerIce Solvers

Github already render .md files but do not support math
You can generate .pdf using e.g. pandoc

To generate the doucumentation of the adjoint method run:
> pandoc -d MakeDoc_Adjoint.yml

## Adjoint methods

 Following a re-organisation of the adjoint methods in Elmer/Ice in April 2020,  

- generic solvers are found under [elmerice/Solvers/Adjoint](../Adjoint)
- model dependent solvers are found under elmerice/Solvers/Adjoint*[Name_Of_Model]*,
where *[Name_Of_Model]* can be Stokes, SSA, Thickness (to come...).

Here is the list of supported and documented solvers:
 
Generic Solvers:  

- [Adjoint_GradientValidation](Adjoint_GradientValidation.md)
- [Adjoint_LinearSolver](Adjoint_LinearSolver.md)
- [Adjoint_CostDiscSolver](Adjoint_CostDiscSolver.md)
- [Adjoint_CostContSolver.md](Adjoint_CostContSolver.md)
- [Adjoint_CostRegSolver.md](Adjoint_CostRegSolver.md)

SSA Solvers:

- [AdjointSSA_SSASolver](AdjointSSA_SSASolver.md)
- [AdjointSSA_GradientSolver](AdjointSSA_GradientSolver.md)
- [AdjointSSA_CostFluxDivSolver](AdjointSSA_CostFluxDivSolver.md)
- [AdjointSSA_CostTaubSolver](AdjointSSA_CostTaubSolver.md)


