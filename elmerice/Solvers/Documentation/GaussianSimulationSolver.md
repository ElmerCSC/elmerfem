## Gaussian Simulation Solver {#Gaussian_simulation}

### General Information
- **Solver Fortran File:** GaussianSimulationSolver.F90
- **Solver Name:** GaussianSimulationSolver
- **Required Input Variable(s):**
  - The *mean* of the distribution as a nodal variable
- **Required Output Variable(s):**
  - The main solver variable is a random sample drawn from the given normal distribution
- **Required Input Keywords:**
  - **Solver Section**:
    - Background Variable name = *String* : Name of the variable that contains the *mean*
    - standard deviation = *Real* : the standard deviation $\sigma$
    - Covariance type = *String*  : Available choices to construct the covariance matrix
      - "diffusion operator"
      - "full matrix"
      - "diagonal"
    - **covariance type specific keywords**: see [CovarianceUtils](#Covariance_Module)
- **Optional Input Keywords:**
    - **Solver Section**:
      - Random Seed = *Integer*: a seed to initialize the random generator for repeatability

### Remark
This documentation contains equations and is part of a generic documentation that can be converted to pdf using pandoc:
```
> pandoc -d MakeDoc_CovarianceUtils.yml
```

### General Description

For a random variable $X$ that is normally distributed as $$X \sim \mathcal{N}(x^b,C)$$, with $x^b$ the mean and $C$ the covariance matrix, it is possible to draw non-conditionnal realizations $x^s$ as
$$x^s = x^b + V.z$$
where:   

- $V$ is obtained from a factorization of $C$ as $C=VV^T$, classically a Cholesky factorisation.  
- $z$ is a vector of uniformly distributed random numbers with zero mean.  

See e.g. Graham et al., (2017).

For an application to uncertainty quantification in ice sheet modeling, using the *diffusion operator* covariance type, see Bulthuis and Larour (2022).

### Implementation

See the generic documentation for [CovarianceUtils](#Covariance_Module) for details on the possible choices to construct the covariance matrix $C$ and for the factorization.

It the solver variable is a vector, each component contains a different realization, otherwise each call to the solver (e.g. during steady-state iterations) will give a different realization.  

### Known Bugs and Limitations

- Limited to serial if using the "full matrix" covariance method.
- The *diffusion operator* might be inaccurate near the boundaries or for highly distorted elements (see Guillet et al., 2019)
- For the moment the implementation is limited to isotropic covariances with spatially uniform parameters (standard deviation and correlation length scale); but this could be improved (see Guillet et al., 2019)

### SIF Contents

```
Solver 1
  Equation = "GSim"
  Variable = -dofs 1 "xs"
  procedure = "ElmerIceSolvers" "GaussianSimulationSolver"

  !# Variable names
  Background Variable Name = String "xb"

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

- Bulthuis, K. and Larour, E.: Implementation of a Gaussian Markov random field sampler for forward uncertainty quantification in the Ice-sheet and Sea-level System Model v4.19, Geosci. Model Dev., 15, 1195–1217, https://doi.org/10.5194/gmd-15-1195-2022, 2022

- Graham, F. S., Roberts, J. L., Galton-Fenzi, B. K., Young, D., Blankenship, D., and Siegert, M. J.: A high-resolution synthetic bed elevation grid of the Antarctic continent, Earth Syst. Sci. Data, 9, 267–279, https://doi.org/10.5194/essd-9-267-2017, 2017.

- Guillet O., Weaver A.T., Vasseur X., Michel Y., Gratton S., Gurol S. Modelling spatially correlated observation errors in variational data assimilation using a diffusion operator on an unstructured mesh. Q. J. R. Meteorol. Soc., 2019. https://doi.org/10.1002/qj.3537
