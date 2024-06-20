## Covariance Vector Multiply Solver {#Covariance_Vector_product}

### General Information
- **Solver Fortran File:** CovarianceVectorMultiplySolver.F90
- **Solver Name:** CovarianceVectorMultiplySolver
- **Required Input Variable(s):**
  - A nodal input variable $x$
- **Required Output Variable(s):**
  - The product of a covariance matrix with the input variable
- **Required Input Keywords:**
  - **Solver Section**:
    - Input Variable  = *String* : Name of the input variable
    - standard deviation = *Real* : the standard deviation $\sigma$
    - Covariance type = *String*  : Available choices to construct the covariance matrix
      - "diffusion operator"
      - "full matrix"
      - "diagonal"
    - **covariance type specific keywords**: see [CovarianceUtils.md](CovarianceUtils.md)
- **Optional Input Keywords:**
    - **Solver Section**:
      - Normalize = *Logical*: wether to normalize the output (default: False)


### General Description

Compute the product $$y=C.x$$
with:   
- $x$ and input variable   
- $C$ a covariance matrix   

Applications:  
- **covariance visualization and code validation**: the spatial correlation function at a given node $z_i$ corresponds to the $i$-th column of the covariance matrix $C$. It can be visualized by plotting the result of applying $C$ to a vector that has a value of one at $z_i$ and a value of zero at all other points (Guillet et al., 2019).  
- **Filtering**: When *Normalize = Logical True*, the output is normalized by the results of applying $C$ to a vector full of ones. If the kernel is a Gaussian correlation function, this would be equivalent to applying a Gaussian filter and this will thus smooth the input variable. The Matérn covariance, obtained with the *diffusion operator* method, converges to the Gaussian correlation function when the smoothness parameters $\nu$ tends to infinity.

### Implementation

See the generic documentation for [CovarianceUtils.md](CovarianceUtils.md) for details on the possible choices to construct the covariance matrix $C$.


### Known Bugs and Limitations

- Limited to 2D meshes.   
- Limited to serial if using the "full matrix" covariance method.   


### SIF Contents

```
Solver 1
 Equation = "Filter"
 Variable = "y"
 Procedure = "CovarianceVectorMultiplySolver" "CovarianceVectorMultiplySolver"

 input variable = string "x"

 Normalize = Logical True !#(default: False)

!# Covariance types
 !############################################################################
 !# keywords for the "diffusion operator" method
 !# see CovarianceUtils.md for other choices
 !############################################################################
  Covariance type = String "diffusion operator"

  Matern exponent m = Integer $m
  correlation range = Real $range
  standard deviation = Real $std

!# The diffusion operator method requires to solve symmetric positive definite
!# linear systems
  Linear System Solver = Direct
  Linear System Direct Method = umfpack

  Linear System Refactorize = Logical False
  Linear System Symmetric = Logical True
  Linear System Positive Definite = Logical True

end
```

### Examples


### References

- Guillet O., Weaver A.T., Vasseur X., Michel Y., Gratton S., Gurol S. Modelling spatially correlated observation errors in variational data assimilation using a diffusion operator on an unstructured mesh. Q. J. R. Meteorol. Soc., 2019. https://doi.org/10.1002/qj.3537
