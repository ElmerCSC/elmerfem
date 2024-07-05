---
title: |
  | ElmerIce Documentation :
  | Spatial correlation modelisation
author:
- F. Gillet-Chaulet
date:
- 20/03/2024
---

## Introduction

This documentation provide the documentation for the following solvers:   
- [BackgroundErrorCostSolver](#Background_Error)   
- [CovarianceVectorMultiplySolver](#Covariance_Vector_product)   
- [GaussianSimulationSolver](#Gaussian_simulation)   

They are all based on operations involving covariances matrices. Different methods have been implemented to compute or apply covariances matrices. These methods are described in the documentation of the module:   
- [CovarianceUtils](#Covariance_Module)   

### Remark
The complete pdf documentation can be obtained using pandoc:   
```
> pandoc -d MakeDoc_CovarianceUtils.yml
```
