# Gradient Validation 

Validate the computation of the gradient of the cost function with the adjoint method against a finite difference 
computation (see [Adjoint_GradientValidation](https://github.com/ElmerCSC/elmerfem/blob/elmerice/elmerice/Solvers/Documentation/Adjoint_GradientValidation.md))

Here we add random noise to the analytical convection velocity and validate the computation of the cost function with respect to the velocity

## Content

* *MakeValidation.sh*: bash script to run several validations 
* *Validation.sif*: .sif input file for the validation

