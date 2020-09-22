## Adjoint GradientValidation

**Module name**: Adjoint_GradientValidation  
**Module subroutines**: Adjoint_GradientValidation  
**Module authors**: Fabien Gillet-Chaulet (IGE-Grenoble)  
**Document authors**: Fabien Gillet-Chaulet  
**Document edited**: 25.04.2020  


### Introduction

This solver is an utility to check the accuracy of the gradients computed from the adjoint solvers.

We have tried to validate each solver independently; however as the model is continuously under developpement,
and the Elmer configuration sometimes complex it might be interesting to check that your derivatives used for
the minimisation are correct.

As the optimisation algorithm uses the gradient for the descent direction, a minimisation that get trapped in linear searches 
or can not decrease the value of the cost function might be the sign for inaccurate gradients.  

> We have done the maximum to derive the adjoints by differentiating crucial parts by hand, however the gradient maight still be 
not as accurate as possible, *e.g.* if the direct problem is non-linear.  
> Also Elmer is continuously improving and allow for a lot
of flexibility, so we can not exclude that some features, *e.g.* different element types, etc..., will not be supported.  
> **Please report if you encounter such problems**.

### Theory
 This solver compare the derivatives computed with the adjoint codes with derivatives obtained by forward finite differences.

The Taylor expansion of the cost function $J(\boldsymbol{x})$ for a small perturbation $h\boldsymbol{x_p}$ leads to
$$ J(\boldsymbol{x} + h\boldsymbol{x_p}) = J(\boldsymbol{x}) + h \boldsymbol{g} . \boldsymbol{x_p} + \mathcal{O}(h||\boldsymbol{x_p}||) $$

This solver then compares the totat derivative $dJ$ computed at a given state $\boldsymbol{x}$ and for a given perturbation 
$\boldsymbol{x_p}$ :
$$ dJ_1=<\boldsymbol{x_p},\boldsymbol{g}> $$
with the forward finite difference equivalent:
$$ dJ_2=\dfrac{J(\boldsymbol{x} + h\boldsymbol{x_p})  - J(\boldsymbol{x})  }{h} $$.

$dJ_1$ is computed at the first iteration, and $dJ_2$ is computed in the following $N$ iterations, starting with $h=1.0$ and
decreasing $h$ by two for each successive iteration.

The relative error is computed as 

$$ \epsilon = \dfrac{abs(dJ_1-dJ_2)}{abs(dJ_1)} $$

If your set-up is correct, $\epsilon$ should tend to $0$.

**Well done!! you can replace the validation solver by the optimisation solver and find your best initial state**

### Keywords

Bellow is the sequence and related keywords in the *.sif* file:  

- First set your number of iterations in the simulation section:
```
Simulation

 Simulation Type = Steady State

 !# set the number of iterations you want to perform
  Steady State Max Iterations = 20
  Steady State Min Iterations = 20

End
```
- Second, set-up your configuration file that computes the cost function and its derivative with respect 
to your *optimisation variable*

- Finally, put the validation solver:

```
Solver id
  Equation = "GradientValidation"
  procedure = "ElmerIceSolvers" "Adjoint_GradientValidation"
 
 !# Name of the variable that contain the cost function value
  Cost Variable Name = String ""
 !# Name for the state variable $x$
  Optimized Variable Name = String ""
 !# Name for the gradient variable $g$
  Gradient Variable Name = String ""
 !# Name of the perturbation variable $xp$
 !# we use $-g$ if not provided.
  Perturbation Variable Name = String ""
 !# Name for the output file
 ! # provide h, relative error, dJ_1, dJ_2
  Result File = File "Validation.dat"

end

```

### Tests and Examples

- See examples for the [Adjoint CostRegSolver](../../examples/Adjoint_CostRegSolver)
