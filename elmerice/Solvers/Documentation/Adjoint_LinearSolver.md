## Adjoint of Linear Systems {#adjoint_linearsolver}

**Module name**: Adjoint_LinearSolver  
**Module subroutines**: Adjoint_LinearSolver  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)  
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 22.04.2020  


### Introduction

This solver is the adjoint code of a solver that solves a linear system and apply dirichlet conditions, 
i.e. it is the adjoint of the following part in a solver:
```fortran
  CALL DefaultDirichletBCs()
  UNorm = DefaultSolve()
```

If the direct problem is non-linear the adjoint should revert the non-linear iterations, which is not the case here;  
However, if the non-linear direct equation is solved using a Newton method, it should be sufficient to solve only the adjoint corresponding to the last linear system solution.

### Theory

When solving a linear system of equation : 

$$ \boldsymbol{A} x = F $$

the adjoint writes: 

> A^T b = xb  

where **x** is the solution of the direct problem, **A^T** is the transpose of the stiffness matrix **A**,
**b** is the *adjoint variable* and **xb** is the *sensitivity* of the cost function with respect to the solution of the direct problem **x**. **xb** should be computed within a cost solver and should be named *velocityb* if the variable of the direct problem is *Flow Solution* (i.e. solving Stokes) or *ssavelocity* (i.e. solving SSA), or *Varb* for other direct equations (where *Var* is the name of the direct solver variable). **xb** has to be computed in the model coordinate system and it will be rotated in Normal-Tangential if needed.

Further the adjoint variables for the stiffness matrix **Ab**, and force vector **Fb** are given by:

> Ab=-bx^T  
> Fb=b

In Elmer, Dirichlet conditions are  applied directly to the matrix structure so that a row is zeroed except for the diagonal which is set to one. Then the r.h.s. value determines the value of the field variable  in the solution of the linear system.  
This will results for the adjoint code to set to 0 the corresping lines in **Ab** and **Fb**, which is equivalent to set **b=0** where Dirichlet conditions have been applied to the direct problem.

The solver will look in the *Body Forces* and *Boundary Conditions* for the keywords associated to Dirichlet conditions for the direct variable **x** and set **b=0** accordingly.

Getting the derivative of the cost function with respect to the direct solver input parameters usually requires to write the adjoint code of the direct solver that fills **A** and **F** and thus depends on the solution **x** and the adjoint variable **b**. This will be solver specific.

### Keywords

Bellow is the sequence and related keywords in the *.sif* file:  

```
!## Direct solver
 !##  compute e.g. the velocity field
  Solver i
   Equation = String "NameOfEquation"
   Variable = String "x" 
   Variable DOFs =  Integer "DOFs"
   ...
  End
  
!## Cost Solver
 !##  compute a cost function, 
 !  e.g. the mismatch between the model velocity and some observations
 !##  compute the sensitivity xb of the cost function 
 !  w.r.t. the solution of the direct solver
 Solver j
   ...
   Should compute the sensitivity "xb" 
   in a variable named:
      "velocityb" (Stokes or SSA) 
      or "Varb" (with Var the name of the solver variable for other solvers)
 End
    
!## Dirichlet conditions for "x" in the Boundary Conditions or Body Forces
 Boundary Condition i
   x_k = ...
   x_k Condition = ....
    
   Normal-Tangential x = Logical True
   Normal-Tangential x condition = ...
 End
   
 Body Force i
   x_k = ...
   x_k Condition = ....
 End 

```
 The solution for the adjoint variable **b** is computed by:

```
 Solver *solver id* 
  
    Equation = String "Adjoint"  
    procedure = "ElmerIceSolvers" "Adjoint_LinearSolver"
 !# The adjoint variable. Can any name as far as it is consistent everywhere 
    Variable = String "Adjoint"  
 !# Degrees of freedom for the adjoint variable, 
 !   needs to be consistent with the DOFs of the Direct solver
    Variable DOFs =  Integer "DOFs"
    
 !# The name of the Equation for the direct Solver; cf above
    Direct Solver Equation Name = String "NameOfEquation"
      
 !# Keywords to solve the linear system
    Linear Sytem Solver = 
      
 End
```
### Limitations and possible improvments
 
- No support for periodic boundary conditions.

### Tests and Examples

- See examples for the [inverse methods](../../examples/Inverse_Methods)
