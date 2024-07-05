## Covariance Utils Module {#Covariance_Module}

### General Information
- **Fortran Module File:** CovarianceUtils.F90

### Remark
This documentation contains equations and is part of a generic documentation that can be converted to pdf using pandoc:
```
> pandoc -d MakeDoc_CovarianceUtils.yml
```

### General Description

The *CovarianceUtils* modules contains routines and functions for the generic operations involving a covariance matrix $B$.
In all the implementation we use the following standard factorisations:
$$B = \Sigma C \Sigma$$
$$B^{-1} = \Sigma^{-1} C^{-1} \Sigma^{-1}$$
where:   

- $C$ is a correlation matrix ,  
- $\Sigma$ is a diagonal standard deviation matrix.

This module is used in the following solvers:   

- [BackgroundErrorCostSolver](#Background_Error)   
- [CovarianceVectorMultiplySolver](#Covariance_Vector_product)   
- [GaussianSimulationSolver](#Gaussian_simulation)   

#### Diffusion operator

The implementation follows Guillet et al. (2019). See Recinos et al. (2023) and Bulthuis and Larour (2022) for similar methods in the context of ice flow modeling.

The motivation of the method is that $C$ is a full-rank matrix of size $(n \times n)$, with n the number of mesh nodes, so that in general for large problem directly building $C$ becomes impossible.

Guillet et al. (2019)  show that applying $m$ successive applications of a discretized representation of the operator $I - l^2\nabla^2$
is equivalent to applying an inverse correlation operator $\mathcal{C}^{-1}$ which as a kernel given by the following Matérn functions:
$$c_{m,l}=\dfrac{2^{2-m}}{(m-2)!} \left ( \dfrac{d}{l} \right )^{m-1} \mathcal{K}_{m-1} (d/l)$$
where:   

- $\mathcal{K}_{m-1}$ is the modified Bessel function of the second kind of order $m$,   
- $d$ is the euclidean distance between 2 points,   
- $l$ is a **correlation lenght scale** (or **range**),   
- $m$ is a **smoothness** parameter that control the *shape* of the function.

Remarks:   

- In the literature, the Matérn functions are often defined with the smoothness parameter $\nu$ which can be a real; Here we are restricted to integers an  $m=\nu+1$.   
- The Matérn functions have two limit cases, the exponential correlation function for $\nu=1/2$ and the squared exponential (or gaussian) correlation function for $\nu \to \infty$.
- The Gaussian limit can be approached by setting the range to the Dayley length scale (Guillet et al. (2019) Eq. 7):   
$$D=\sqrt{2m-4}l$$

For the following, we define the following the mass matrix $M$ and stiffness matrix $K$ discretized by the FEM:
$$M_{ij}=\int_{\Omega} \phi_i \phi_j d\Omega$$
$$K= l^2 \int_{\Omega} \nabla \phi_i \nabla \phi_j d\Omega$$

The [BackgroundErrorCostSolver](#Background_Error) requires the inverse correlation matrix $C^{-1}$ which is discretized as:
$$C^{-1} = \Gamma^{-1}ML^{-1}_M \Gamma^{-1}$$
with:   

- $\Gamma=\sqrt{4\pi(m-1)}l I$ is a normalization matrix      
- $L^{-1}_M = [M^{-1}(M + K)]^m$

The [CovarianceVectorMultiplySolver](#Covariance_Vector_product) requires the correlation matrix $C$ which is discretized as:
$$C = \Gamma L_M M^{-1} \Gamma$$
with:  
 $$L_M = [(M + K)^{-1}M]^m$$

The [GaussianSimulationSolver](#Gaussian_simulation) requires a square root of $B$ which is obtained from the following factorisation:
$$B=VV^T$$
with:   
$$V= \Sigma \Gamma L^{1/2}_M (M^{-1/2})^T$$
where:   

- to compute $M^{-1/2}$ we take the *lumped* mass matrix,   
- $L^{1/2}_M = [(M + K)^{-1}M]^{m/2}$, so restricting the application to even values of $m$

**Remark:** This discretization implies Neumann boundary conditions and it is known that it is inaccurate near the boundaries; cf sec. 3.6 and 5.4 in Guillet et al. (2019).   

Keywords related to this method:   

- the method is chosen with the keyword *Covariance type = String "diffusion operator"*   
- the correlation length scale $l$ is assumed uniform and given by the keyword *correlation range = Real*   
- the smoothness parameter $m$ is given by the keyword *Matern exponent m = Integer* ($m>2$ for all solvers and m must be even for the GaussianSimulationSolver)   
- the standard deviation to build $\Sigma$ is assumed uniform and given by the keyword *standard deviation = Real*

Limitations:   

- can be used in serial/parallel.   
- as only be tested on 2D meshes or on a 2D boundary of a 3D mesh. In the latter case *projection coordinate=Integer* sets the corresponding coordinates to 0, so that the operations are performed in the projected plane; e.g. if *projection coordinate= Integer 3*, the *z-* coordinates are set to 0 to compute the mass and stiffness matrices


sif example:  
```
Solver 1
...
 !############################################################################
 !# Covariance types
 !############################################################################
   Covariance type = String "diffusion operator"

  Matern exponent m = Integer $m
  correlation range = Real $range
  standard deviation = Real $std

  !# Whene used as a boundary solver,
  !#projection coordinate = Integer ...

!# The diffusion operator method requires to solve symmetric positive definite
!# linear systems;
  Linear System Solver = Direct
  Linear System Direct Method = umfpack

  Linear System Refactorize = Logical False
  Linear System Symmetric = Logical True
  Linear System Positive Definite = Logical True

end
```


#### Full matrix

This has been implement mainly for testing/validation on small serial test cases. Here the we build the matrix $C$ using standard analytical correlation functions.

Operations on $C$ are performed using standard Lapack linear algebra routines:   

- For the [GaussianSimulationSolver](#Gaussian_simulation) with use a Cholesky decomposition using the Lapack routine *dpptrf*   
- For the [BackgroundErrorCostSolver](#Background_Error), the inverse is obtained using the lapack routine *dpptri* after the Cholesky factorisation.

The following correlation functions have been implemented:

- the exponential:
$$c(d)=exp(-d/l)$$   
- the Gaussian correlation function:
$$c(d)=exp(-d^2/(2l^2))$$
- the Matérn function for integer values of $\nu$ (*materni*):
$$c(d)=\dfrac{2^{1-\nu}}{(\nu-1)!} \left ( \dfrac{d}{l} \right )^{\nu} \mathcal{K}_{\nu} (d/l)$$
- Analytical power solution of the Matérn function when $\nu$ is half integer so that  $p=\nu-1/2$ is an integer (*maternp*) (rq. p=0 is the exponential):   
  - p=1
  $$c(d)=(1+d/l)exp(-d/l)$$   
  - p=2
  $$c(d)=(1+d/l+d^2/(3l^2))exp(-d/l)$$


**Remark** for *MaternI* the Bessel function $\mathcal{K}$ is obtained by recursion and accuracy for high values of $n$ ($n>10$) has not been tested.

Keywords related to this method:

- the method is chosen with the keyword *Covariance type = String "full matrix"*   
- the correlation length scale $l$ is assumed uniform and given by the keyword *correlation range = Real*   
- The correlation function is given by the keyword *correlation type = String*; possible choices are:    

      - *"exponential"*      
      - *"gaussian"*   
      - *"materni"*   
      - *"maternp"*       

- For *materni* the integer value of $\nu$ is given by *MaternI order = Integer*      
- For *maternp* the integer value of $p$ is given by *MaternP order = Integer* (restricted to 1 and 2)   


Limitations:

- can only be used in serial.   
- restricted to relatively small problems      
- Can be used as a boundary solver; in this case if the mesh is 3D, only the $x$ and $y$ coordinates are used to compute the euclidean distance.   


sif example:
```
Solver 1
...
 !############################################################################
 !# Covariance types
 !############################################################################
   Covariance type = String "full matrix"

  correlation range = Real $range
  standard deviation = Real $std

  correlation type = Sting "exponential"
  #or
  # correlation type = Sting "gaussian"
  #or
  # correlation type = Sting "materni"
  # MaternI order = Integer 1
  #or
  # correlation type = Sting "maternp"
  # MaternI order = Integer 1


  # rq. there is no linear system to solve
end
```

#### Diagonal

This is the simple choice, there is no spatial correlation and the covariance matrix is simply $$B=\sigma^2 I$$

Keywords related to this method:

- the method is chosen with the keyword *Covariance type = String "diagonal"*

```
Solver 1
...
 !############################################################################
 !# Covariance types
 !############################################################################
   Covariance type = String "diagonal"

   standard deviation = Real $std
```

### Module functions and subroutines

Generic routines for the initialisation of the matrices (*D* referring to the *diffusion operator* and *L* to the *full matrix*):
```
INTERFACE CovarianceInit
  MODULE PROCEDURE CovarianceInitD,CovarianceInitL
END INTERFACE
```

Generic routines to perform covariance matrix vector multiplications:   
```
INTERFACE CovarianceVectorMultiply
  MODULE PROCEDURE CovarianceVectorMultiplyD,CovarianceVectorMultiplyL
END INTERFACE
```

Generic routines to perform inverse covariance matrix vector multiplications:   
```
INTERFACE InvCovarianceVectorMultiply
  MODULE PROCEDURE InvCovarianceVectorMultiplyD,InvCovarianceVectorMultiplyL
END INTERFACE
```
Generic routines to perform square-root covariance matrix vector multiplications:   
```
INTERFACE SqrCovarianceVectorMultiply
   MODULE PROCEDURE SqrCovarianceVectorMultiplyD,SqrCovarianceVectorMultiplyL
END INTERFACE
```

Functions related to the computation of the analytical correlation functions.

### Examples

- ElmerIce unitary tests:
   - [ELMER_TRUNK]/elmerice/Tests/CovarianceVector
   - [ELMER_TRUNK]/elmerice/Tests/CovarianceVector2

- Validation test cases and set-ups:   
   - https://gricad-gitlab.univ-grenoble-alpes.fr/gilletcf/CovarianceUtils

### References

- Recinos, B., Goldberg, D., Maddison, J. R., and Todd, J.: A framework for time-dependent ice sheet uncertainty quantification, applied to three West Antarctic ice streams, The Cryosphere, 17, https://doi.org/10.5194/tc-17-4241-2023, 2023.
- Bulthuis, K. and Larour, E.: Implementation of a Gaussian Markov random field sampler for forward uncertainty quantification in the Ice-sheet and Sea-level System Model v4.19, Geosci. Model Dev., 15, 1195–1217, https://doi.org/10.5194/gmd-15-1195-2022, 2022
- Guillet O., Weaver A.T., Vasseur X., Michel Y., Gratton S., Gurol S. Modelling spatially correlated observation errors in variational data assimilation using a diffusion operator on an unstructured mesh. Q. J. R. Meteorol. Soc., 2019. https://doi.org/10.1002/qj.3537
