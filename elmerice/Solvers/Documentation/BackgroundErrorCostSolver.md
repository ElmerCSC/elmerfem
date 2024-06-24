## Background Error Cost Solver {#Background_Error}

### General Information
- **Solver Fortran File:** BackgroundErrorCostSolver.F90
- **Solver Name:** BackgroundErrorCostSolver
- **Required Input Variable(s):**
  - Variable: the *active* variable
  - Background Variable: A *prior* estimate against which we compare the *active* variable
- **Required Output Variable(s):**
  - Cost Variable: global variable containing the *cost*
  - Gradient Variable: derivative of the cost w.r.t. the *active* variable
- **Required Input Keywords:**
  - **Solver Section**:
    - Variable name = *String* : Name of the *active* variable
    - Background Variable name  = *String* : Name of the *Background Variable*
    - Gradient Variable name = *String* : Name of the *Gradient Variable*
    - Cost Variable name  = *String* : Name of the *Cost Variable*
    - Cost Filename = *File* : Cost file name
    - Reset Cost Value = *Logical* (Default: True) : Cost and gradient are initialized to 0 if true
    - standard deviation = *Real* : the standard deviation $\sigma$
    - Covariance type = *String*  : Available choices to construct the covariance matrix
      - "diffusion operator"
      - "full matrix"
      - "diagonal"
    - **covariance type specific keywords**: see [CovarianceUtils](#Covariance_Module)

### Remark
This documentation contains equations and is part of a generic documentation that can be converted to pdf using pandoc:
```
> pandoc -d MakeDoc_CovarianceUtils.yml
```

### General Description

This solver is mainly intended to be used with the adjoint inverse methods implemented in Elmer/Ice (see the corresponding documentations).

It is an alternative to the *Regularisation* solver that penalizes the fist spatial derivatives of the *active* variable $x$, where $x$ is usually the basal friction coefficient, ice viscosity or bed elevation (see [Adjoint_CostRegSolver.md](Adjoint_CostRegSolver.md))

Here the cost is computed as $$Cost=0.5 (x-x^b). B^{-1} .(x-x^b)$$ where:

- $x$ is the *active* variable   
- $x^b$ is the *background* variable, i.e. a prior estimate of $x$   
- $B$ is the background error **covariance** matrix, i.e. is described the statics of the *expected* (or *tolerated*) difference between the background $x^b$ and the active variable $x$.

The gradient of the Cost w.r.t. the *active* variable is then obtained as:
$$dCost/dx=B^{-1}.(x-x_b)$$


This cost function is usually applied in addition of a cost function that penalizes the difference between a diagnostic variable $u(x)$ and its observed counterpart $u^o$, which is general would be written as $$Cost^o=0.5 (u(x)-u^o). R^{-1} .(u(x)-u^o)$$, where $R$ is the observation error correlation matrix. This cost can be computed using e.g. the [Adjoint_CostDiscSolver](Adjoint_CostDiscSolver.md) (restricted to errors that are not spatially correlated, i.e. $R$ is diagonal and contains the observation errors variances)

If the observation and background errors are gaussian (described by the respective covariance matrices $R$ and $B$), the errors unbiaised and the observation operator is linear, i.e. $u=K.x$, the solution of the variational inverse problem, i.e. finding the minimum of the cost function, solves the Bayesian inference problem.  The optimal solution then depends on the confidence that is given to the observation and to the background, which is *parameterized* by the corresponding covariance matrices. Special care must then be taken when prescribing these matrices.

Providing:   

- the input active variable $x$, given by the solver keyword *Variable name*  
- a variable containing the background $x^b$, given by the solver keyword *Background Variable name*   
- a method and the parameters for the covariance structure, given by the solver keyword *Covariance type* and method specific keywords   
this solver computes:   
- the *Cost* with is saved in the global cost variable, given by the solver keyword *Cost Variable name*      
- the gradient of the cost variable w.r.t. $x$, given by the solver keyword *Cost Variable name*     

The value of the cost a a function of the iteration is save in an ascii file defined by the solver keyword *Cost Filename*.

### Implementation

The covariance matrix is of size $n \times n$, with $n$ the number of mesh nodes, is usually full rank, and often poorly known. It is then standard to parameterize this matrix using standard correlation functions $c(d)$ that describe the spatial correlation between 2 points as a function of the distance $d$ between the points. The correlation structure depends then on the correlation function (typical functions are, e.g., the exponential, squared exponential (or gaussian), Matérn functions, see [CovarianceUtils](#Covariance_Module)) which usually depends on a *correlation length scale* or *range*.

The inverse covariance matrix $B^{-1}$ is factored in the standard form
$$B^{-1}=\Sigma^{-1}. C^{-1} . \Sigma^{-1}$$
where:   

- $\Sigma$ is a diagonal matrix containing the  standard deviations, assumed spatially uniforms fro now, i.e. $\Sigma=\sigma I$, with $I$ the identity matrix.   
- $C^{-1}$ is the inverse of the  correlation matrix whose components are defined using standard correlation functions $c(d)$.

In general it is not necessary to explicitly compute and store $B$ (or its inverse), and it can be replaced by an equivalent operator. For  **Covariance type = String "diffusion operator"**, the operator results from the discretization of a diffusion equation, which can be done efficiently for unstructured meshes with the finite element method.  With this method the operator kernel is a correlation function from the Matérn family. The implementation follows  Guillet et al. (2019).

See [CovarianceUtils](#Covariance_Module) for details on the possible choices to construct $C$.

### Discussion

 Brasseur et al. (1996) have shown that adding a smoothness constraint that penalizes a combination of the nom and of the spatial derivatives up to order 2, is equivalent, for an infinite domain, to imposing a kernel from the  Matérn family with a **smoothness parameter** $\nu=1$. This has been generalized to higher dimensions and derivatives by Barth et al. (2014). Regularisation of inverse problems can often be reinterpreted in the Bayesian framework (Calvetti and Somersalo, 2018), so that the effect of this solver will be similar to the classically used *Regularisation* solver that penalizes he first spatial derivatives, and the choice of the correlation structure and parameters will control the **smoothness** of the inverted field. However this solver is then much more versatile and the parameters have a direct physical interpretation.

 For an application of this method in ice-sheet modeling for the inversion of both basal friction and viscosity in the Antarctic Ice Sheet see e.g. Recinos et al. (2023).


### Known Bugs and Limitations

- Limited to serial if using the "full matrix" covariance method.   
- The *diffusion operator* might be inaccurate near the boundaries or for highly distorted elements (see Guillet et al., 2019)
- For the moment the implementation is limited to isotropic covariances with spatially uniform parameters (standard deviation and correlation length scale); but this could be improved (see Guillet et al., 2019)

### SIF Contents


```
Solver 1
  Equation = String "CostReg"
  procedure = "ElmerIceSolvers" "BackgroundErrorCostSolver"
  Variable = -nooutput "dumy"

  !# Variable names
  Variable Name = String "bed"
  Gradient Variable Name = String "bedb"
  Background Variable Name = String "bmean"
  Cost Variable Name= String "CostValue"

  !# output cost file
  Cost Filename = File "CostFile.dat"

  !# True if cost function and gradient must be initialized to 0 in this solve
  !# otherwise cost and gradient will be added to the previous Values
  !# which is the case if this solver is used after a Cost Solver
  !# measuring error w.r.t. observations
  Reset Cost Value = Logical False

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

- ElmerIce unitary tests:     
   - [ELMER_TRUNK]/elmerice/Tests/BackgroundError

### References

- Barth, A., Beckers, J.-M., Troupin, C., Alvera-Azcárate, A., and Vandenbulcke, L.: divand-1.0: n-dimensional variational data analysis for ocean observations, Geosci. Model Dev., 7, 225–241, https://doi.org/10.5194/gmd-7-225-2014, 2014.
- Brasseur, P., Beckers, J.M.,   Brankart, J.M. and R. Schoenauen, Seasonal temperature and salinity fields in the Mediterranean Sea: Climatological analyses of a historical data set, Deep Sea Research Part I: Oceanographic Research Papers, 43(2), 1996, https://doi.org/10.1016/0967-0637(96)00012-X
- Calvetti D, Somersalo E. Inverse problems: From regularization to Bayesian inference. WIREs Comput Stat., 2018 https://doi.org/10.1002/wics.1427
- Guillet O., Weaver A.T., Vasseur X., Michel Y., Gratton S., Gurol S. Modelling spatially correlated observation errors in variational data assimilation using a diffusion operator on an unstructured mesh. Q. J. R. Meteorol. Soc., 2019. https://doi.org/10.1002/qj.3537
- Recinos, B., Goldberg, D., Maddison, J. R., and Todd, J.: A framework for time-dependent ice sheet uncertainty quantification, applied to three West Antarctic ice streams, The Cryosphere, 17, https://doi.org/10.5194/tc-17-4241-2023, 2023.
